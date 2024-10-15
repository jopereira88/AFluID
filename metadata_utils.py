#!/usr/bin/python3
from abc import ABC, abstractmethod
import re
import pickle


class MetadataTable(ABC):
    '''
    Abstract class to handle the project's metadata tables
    '''
    def __init__(self, filename:str,delim:str, key_index:int = 0):
        '''
        Class initializer
        Accepts: filename (str) - path to te table
        delim (str) - the tabular delimitator
        '''
        self.filename=filename
        self.delim=delim
        self.key_index=key_index
        self.data = self.load_data()
        self.headers = self.load_headers()
    
    def load_data(self):
        '''
        Loads data from the tabular file and stores the data 
        in a dictionary of lists
        '''
        data={}
        try:
            with open(self.filename,'r') as table:
                lines=table.readlines()
                for i in range(len(lines)):
                    if i != 0:
                        line=lines[i].strip()
                        line=line.split(self.delim)
                        data[line[self.key_index]]=[]
                        for i in range(len(line)):
                            if i != self.key_index:
                                data[line[self.key_index]].append(line[i])

            print(f'Data loaded from {self.filename}, successfully.')
        except Exception as e:
            print(f'Failed to load data from {self.filename}: {e}')
            raise 
        return data
    
    def load_headers(self):
        '''
        Loads headers from the tabular file and stores the data 
        in a dictionary
        '''
        try:
            with open(self.filename,'r') as table:
                lines=table.readlines()
                headers={name: index for index, name in enumerate(lines[0].strip().upper().split(self.delim)) if index < 14}

        except Exception as e:
            print(f'Failed to load data from {self.filename}: {e}')
            raise e
        return headers

    @abstractmethod
    def process_metadata(self):
        '''
        Abstract method to process metadata
        Must be implemmented by subclasses
        '''
        pass

    @abstractmethod
    def validate_metadata(self):
        '''
        Abstract method to validate metadata
        Must be implemmented by subclass
        '''
        pass

    @abstractmethod
    def update_metadata(self):
        '''
        Abstract method to update metadata
        Must be implemmented by subclass
        '''
        pass

    @abstractmethod
    def export_metadata(self):
        '''
        Abstract method to export metadata
        Must be implemmented by subclass
        '''
        pass


class SequenceMetadata(MetadataTable):
    def __init__(self, filename: str):
        super().__init__(filename, delim=',', key_index=0)
    
    def validate_metadata(self):
        '''
        Returns the number of empty fields on each column as a dictionary
        '''
        empty={col : 0 for col in self.headers.keys() if self.headers[col] != self.key_index }
        for row in self.data.values():
            for col in self.headers:
                if self.headers[col] != self.key_index and row[self.headers[col]-1]=='':
                    empty[col]+=1
        return empty
    
    def process_metadata(self):
        return super().process_metadata()
    
    def update_metadata(self,key:str,field:str,value:str):
        '''
        Updates metadata values using key list and column to replace values
        '''
        index=self.headers[field.upper()] - 1
        for i in self.data:
            if i == key:
                self.data[key][index]=value
    
    def get_headers(self):
        '''
        Returns the name of the headers
        '''
        headers=[]
        for key in self.headers:
            if self.headers[key]==self.key_index:
                headers.append(f'KEY:{key}')
            else:
                headers.append(key)
        return headers

    def get_values(self, field:str):
        '''
        Get a value per row of a certain field in self.headers
        Returns a dictionary of not null values and the keys
        '''
        index=self.headers[field.upper()] - 1
        query={}
        for key in self.data:
            if self.data[key][index] != '':
                query[key]=self.data[key][index]
        return query
    
    def get_keys_per_notnull(self,field:str):
        '''
        Returns a list of keys that have a certain not null value
        as field.
        '''
        keys=[key for key in self.data if self.data[key][self.headers[field.upper()]-1] != '']
        return keys
    
    def get_keys_per_2fields(self,field1:str,field2:str,value1:str,value2:str=''):
        '''
        returns a list of keys that have a not null value1 on field1 
        and a set value on field2 
        '''
        keys=[key for key in self.data if (self.data[key][self.headers[field1.upper()]-1] == value1 and
              self.data[key][self.headers[field2.upper()]-1] == value2)]
        return keys
    
    def export_metadata(self,filename:str):
        '''
        exports metadata to a .csv file tith a specified filename
        and filepath
        '''
        with open(filename,'w') as table:
            for header in self.headers.keys():
                table.write(f'{header}{self.delim}')
            table.write('\n')
            for key in self.data:
                table.write(f'{key}{self.delim}')
                for i in range(len(self.data[key])):
                    table.write(f'{self.data[key][i]}{self.delim}')
                table.write('\n')

class ClusterMetadata(MetadataTable):
    def __init__(self, filename: str):
        super().__init__(filename, delim="\t" , key_index=0)
        self.process_metadata()
    
    def validate_metadata(self ,field:str):
        '''
        Returns the number of not singular value clusters of a 
        declared field with multiple values
        '''
        index=self.headers[field.upper()] - 1
        multi_value={}
        if index > 0 and index < len(self.headers)-1:
            for key in self.data:
                if len(self.data[key][index]) >1:
                    multi_value[key]=self.data[key][index]

        return multi_value
    
    def process_metadata(self):
        '''
        Transforms the multi-key fields in dictionaries  
        '''
        for key in self.data:
             for i in range(1,len(self.data[key])):
                self.data[key][i]=eval(self.data[key][i])
    
    def update_metadata(self):
        return super().update_metadata()
    
    def get_headers(self):
        '''
        Returns the name of the headers
        '''
        headers=[]
        for key in self.headers:
            if self.headers[key]==self.key_index:
                headers.append(f'KEY:{key}')
            else:
                headers.append(key)
        return headers
    def get_field(self,field:str):
        '''
        Returns a dict with key and field values
        '''
        index=self.headers[field.upper()] - 1
        col={key:value for key in self.data for value in self.data[key][index]}
        return col

    def get_H(self):
        '''
        Returns H genotype from segment 4 clustes 
        that have no mixed genotypes and prints to stdout the mixed gen clusters
        '''
        clust_gen4={}
        mixed=[]
        for key in self.data:
            for seg in self.data[key][self.headers['segments'.upper()]-1].keys():
                if seg=='4':
                    for gen in self.data[key][self.headers['genotypes'.upper()]-1].keys():
                        if gen == 'mixed':
                            mixed.append(key)
                        else:
                            patt=re.compile(r'H\d+',re.IGNORECASE)
                            ex=re.search(patt, gen)
                            if ex is not None:
                                clust_gen4[key]=set()
                                clust_gen4[key].add(ex.group())
        with open('mixedHclusters.pkl','wb') as pkl:
            pickle.dump(mixed,pkl)
        return clust_gen4
    
    def get_N(self):
        '''
        Returns H genotype from segment 4 clustes 
        that have no mixed genotypes and prints to stdout the mixed gen clusters
        '''
        clust_gen6={}
        mixed=[]
        for key in self.data:
            for seg in self.data[key][self.headers['segments'.upper()]-1].keys():
                if seg=='6':
                    for gen in self.data[key][self.headers['genotypes'.upper()]-1].keys():
                        if gen == 'mixed':
                            mixed.append(key)
                        else:
                            patt=re.compile(r'N\d+',re.IGNORECASE)
                            ex=re.search(patt, gen)
                            if ex is not None:
                                clust_gen6[key]=set()
                                clust_gen6[key].add(ex.group())
        with open('mixedNclusters.pkl','wb') as pkl:
            pickle.dump(mixed,pkl)
        return clust_gen6
    def query_clust(self,field:str,query:list):
        '''
        Outputs field information from cluster query list
        '''
        index=self.headers[field.upper()] - 1
        for item in query:
            print(item,self.data[item][index])
    def convert_to_prop(self,field):
        '''
        converts a multi-value field of absolute counts to proprotions
        works on fields 1-3
        '''
        index=self.headers[field.upper()] - 1
        perc={}
        for key in self.data:
            elems=[]
            elems.append(self.data[key][index])
            tot=0
            for dic in elems:
                tot=sum(list(dic.values()))
                for i in dic:
                    dic[i]=round(dic[i]/tot,3)
                perc[key]=dic
        return perc
    
    def export_metadata(self):
        return super().export_metadata()

class BlastReportTable(MetadataTable):
    def __init__(self, filename:str):
        super().__init__(filename, delim='\t',key_index=0)

    def process_metadata(self):
        return super().process_metadata()
    
    def validate_metadata(self):
        return super().validate_metadata()
    
    def update_metadata(self):
        return super().update_metadata()
    
    def get_segment(self):
        '''
        Returns segment for each query acession number in a dict
        '''
        id_seg={}
        field='segment'
        index=self.headers[field.upper()] - 1
        try:    
            for key in self.data:
                id_seg[key]=self.data[key][index]
        except IndexError:
            id_seg[key]='NA'
        return id_seg
    def get_genotype(self):
        '''
        Returns segment for each query acession number in a dict
        '''
        id_gen={}
        field='genotype'
        index=self.headers[field.upper()] - 1
        try:    
            for key in self.data:
                id_gen[key]=self.data[key][index]
            
        except IndexError:
            id_gen[key]='NA'
        return id_gen
    def get_host(self):
        '''
        Returns segment for each query acession number in a dict
        '''
        id_host={}
        field='host'
        index=self.headers[field.upper()] - 1
        try:    
            for key in self.data:
                id_host[key]=self.data[key][index]
        except IndexError:
            id_host[key]='NA'
        return id_host
    def export_metadata(self):
        return super().export_metadata()
    
    def extract_h(self):
        '''
        Extracts the Hx genotype from the genotype field and returns a
        dict with a set of Hx or Hx+mixed 
        '''
        to_mine=self.get_genotype()
        h_reg=re.compile(r'H\d+',re.IGNORECASE)
        mixed_reg=re.compile(r'mixed',re.IGNORECASE)
        h_dict={}
        for key in to_mine:
            h_dict[key]=set()
            if type(to_mine[key])==str:
                ex=re.search(h_reg,to_mine[key])
                ex2=re.search(mixed_reg,to_mine[key])
                
                if ex is not None:
                    h_dict[key].add(ex.group())
                if ex2 is not None:
                    h_dict[key].add(ex2.group())
            else:
                for value in to_mine[key]:
                    ex=re.search(h_reg,value)
                    ex2=re.search(mixed_reg,value)
                    
                    if ex is not None:
                        h_dict[key].add(ex.group())
                    if ex2 is not None:
                        h_dict[key].add(ex2.group())
        return h_dict 
    def extract_n(self):
        '''
        Extracts the Hx genotype from the genotype field and returns a
        dict with a set of Hx or Hx+mixed 
        '''
        to_mine=self.get_genotype()
        n_reg=re.compile(r'N\d+',re.IGNORECASE)
        mixed_reg=re.compile(r'mixed',re.IGNORECASE)
        n_dict={}
        for key in to_mine:
            n_dict[key]=set()
            if type(to_mine[key])==str:
                ex=re.search(n_reg,to_mine[key])
                ex2=re.search(mixed_reg,to_mine[key])
                
                if ex is not None:
                    n_dict[key].add(ex.group())
                if ex2 is not None:
                    n_dict[key].add(ex2.group())
            else:
                for value in to_mine[key]:
                    ex=re.search(n_reg,value)
                    ex2=re.search(mixed_reg,value)
                    
                    if ex is not None:
                        n_dict[key].add(ex.group())
                    if ex2 is not None:
                        n_dict[key].add(ex2.group())
        return n_dict 

class ClusterReportTable(MetadataTable):
    def __init__(self, filename:str):
        super().__init__(filename, delim='\t',key_index=0)
        self.process_metadata()

    def process_metadata(self):
        '''
        Transforms the multi-key fields in dictionaries  
        '''
        for key in self.data:
             for i in range(3,len(self.data[key])):
                self.data[key][i]=eval(self.data[key][i])

    def validate_metadata(self):
        return super().validate_metadata()
    
    def update_metadata(self):
        return super().update_metadata()
    
    def get_segment(self):
        '''
        Returns segment value for each query acession number in a dict
        '''
        id_seg={}
        
        field='segment'
        index=self.headers[field.upper()] - 1
        for i in self.data:
            if self.data[i][0] !='NA':
                keys=[]
                for seg in self.data[i][index].keys():
                    keys.append(seg)
                if len(keys)>1:
                    id_seg[i]=keys
                else:
                    id_seg[i]=keys[-1] 


        return id_seg
    def get_genotype(self):
        '''
        Returns segment for each query acession number in a dict
        '''
        id_gen={}
        field='genotypes'
        index=self.headers[field.upper()] - 1
        for i in self.data:
            if self.data[i][0] !='NA':
                keys=[]
                for seg in self.data[i][index].keys():
                    keys.append(seg)
                if len(keys)>1:
                    id_gen[i]=keys
                else:
                    id_gen[i]=keys[-1] 

        return id_gen
    def get_host(self):
        '''
        Returns segment for each query acession number in a dict
        '''
        id_host={}
        field='hosts'
        index=self.headers[field.upper()] - 1
        for key in self.data:
            if self.data[key][0]!='NA':
                id_host[key]=self.data[key][index]
        
        return id_host
    def print_report(self):
        '''
        Prints report to stdout
        '''
        unassigned=[]
        for key in self.headers:
            print(key,end='\t')
        print()
        for key in self.data:
            if self.data[key][0]!='NA':
                print(key, end='\t')
                for i in range(len(self.data[key])):
                    print(self.data[key][i],end='\t')
                print()
            else:
                unassigned.append(key)
        print('Unassigned accessions:')
        for i in unassigned:
            print(i)

    def extract_h(self):
        '''
        Extracts the Hx genotype from the genotype field and returns a
        dict with a set of Hx or Hx+mixed 
        '''
        to_mine=self.get_genotype()
        h_reg=re.compile(r'H\d+',re.IGNORECASE)
        mixed_reg=re.compile(r'mixed',re.IGNORECASE)
        h_dict={}
        for key in to_mine:
            h_dict[key]=set()
            if type(to_mine[key])==str:
                ex=re.search(h_reg,to_mine[key])
                ex2=re.search(mixed_reg,to_mine[key])
                
                if ex is not None:
                    h_dict[key].add(ex.group())
                if ex2 is not None:
                    h_dict[key].add(ex2.group())
            else:
                for value in to_mine[key]:
                    ex=re.search(h_reg,value)
                    ex2=re.search(mixed_reg,value)
                    
                    if ex is not None:
                        h_dict[key].add(ex.group())
                    if ex2 is not None:
                        h_dict[key].add(ex2.group())
        return h_dict 
    def extract_n(self):
        '''
        Extracts the Nx genotype from the genotype field and returns a
        dict with a set of Nx or Nx+mixed 
        '''
        to_mine=self.get_genotype()
        n_reg=re.compile(r'N\d+',re.IGNORECASE)
        mixed_reg=re.compile(r'mixed',re.IGNORECASE)
        n_dict={}
        for key in to_mine:
            n_dict[key]=set()
            if type(to_mine[key])==str:
                ex=re.search(n_reg,to_mine[key])
                ex2=re.search(mixed_reg,to_mine[key])
                
                if ex is not None:
                    n_dict[key].add(ex.group())
                if ex2 is not None:
                    n_dict[key].add(ex2.group())
            else:
                for value in to_mine[key]:
                    ex=re.search(n_reg,value)
                    ex2=re.search(mixed_reg,value)
                    
                    if ex is not None:
                        n_dict[key].add(ex.group())
                    if ex2 is not None:
                        n_dict[key].add(ex2.group())
        return n_dict 
    
    def convert_to_prop(self,field):
        '''
        converts a multi-value field of absolute counts to proprotions
        works on fields 4-6
        '''
        index=self.headers[field.upper()] - 1
        perc={}
        for key in self.data:
            try:
                elems=[]
                elems.append(self.data[key][index])
                tot=0
                for dic in elems:
                    tot=sum(list(dic.values()))
                    for i in dic:
                        dic[i]=f'{round(dic[i]/tot,5)*100}%'
                    perc[key]=dic
            except IndexError:
                continue
        return perc


    def export_metadata(self):
        return super().export_metadata()


##############################################################################################################################################


""" clust_rep=ClusterReportTable('reports/rand_100_0_report_compiled.txt')
clust_rep.convert_to_prop('genotypes')
for key in clust_rep.data:
    if clust_rep.data[key][0] !='NA':
        print(clust_rep.data[key][3]) """
'''
hemagglutinin=metadata.get_H()
neuraminidase=metadata.get_N()
with open('mixedHclusters.pkl','rb') as pkl:
    H_mixed=pickle.load(pkl)
with open('mixedNclusters.pkl','rb') as pkl:
    N_mixed=pickle.load(pkl)
puremixed={'mixed':0,'assigned':0}
mixed=[]
for cluster in H_mixed:
    if cluster in hemagglutinin.keys():
        puremixed['assigned']+=1
        mixed.append(cluster)
    else:
        puremixed['mixed']+=1
'''