#!/usr/bin/python3
from abc import ABC, abstractmethod
import re
import json
import ast
from collections import namedtuple
from typing import Any


class MetadataTable(ABC):
    """Abstract base class for project metadata tables."""
    def __init__(self, filename:str,delim:str, key_index:int = 0) -> None:
        """
        Initialize a metadata table wrapper.

        Parameters
        ----------
        filename : str
            Path to the source table.
        delim : str
            Field delimiter used by the table.
        key_index : int, default 0
            Index of the column used as the primary key.
        """
        self.filename=filename
        self.delim=delim
        self.key_index=key_index
        self.data = self.load_data()
        self.headers = self.load_headers()
    
    def load_data(self) -> dict[str, list[str]]:
        """
        Load rows from the source table.

        Returns
        -------
        dict
            Mapping from row key to the remaining field values.
        """
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

            #print(f'Data loaded from {self.filename}, successfully.')
        except Exception as e:
            print(f'Failed to load data from {self.filename}: {e}')
            raise 
        return data
    
    def load_headers(self) -> dict[str, int]:
        """
        Load header names from the source table.

        Returns
        -------
        dict
            Mapping from uppercase header names to column indexes.
        """
        try:
            with open(self.filename,'r') as table:
                lines=table.readlines()
                headers={name: index for index, name in enumerate(lines[0].strip().upper().split(self.delim)) if index < 14}

        except Exception as e:
            print(f'Failed to load data from {self.filename}: {e}')
            raise e
        return headers

    @abstractmethod
    def process_metadata(self) -> None:
        """Process loaded metadata in place."""
        pass

    @abstractmethod
    def validate_metadata(self) -> None:
        """Validate loaded metadata."""
        pass

    @abstractmethod
    def update_metadata(self) -> None:
        """Update metadata values."""
        pass

    @abstractmethod
    def export_metadata(self) -> None:
        """Export metadata to a persistent representation."""
        pass


#########################################################################################################################


class SequenceMetadata(MetadataTable):
    def __init__(self, filename: str) -> None:
        """
        Initialize a sequence metadata table backed by a CSV file.

        Parameters
        ----------
        filename : str
            Path to the sequence metadata CSV file.
        """
        super().__init__(filename, delim=',', key_index=0)
    
    def validate_metadata(self) -> dict[str, int]:
        """
        Count empty values per metadata column.

        Returns
        -------
        dict
            Mapping from column name to the number of empty cells.
        """
        empty={col : 0 for col in self.headers.keys() if self.headers[col] != self.key_index }
        for row in self.data.values():
            for col in self.headers:
                if self.headers[col] != self.key_index and row[self.headers[col]-1]=='':
                    empty[col]+=1
        return empty
    
    def process_metadata(self) -> None:
        """
        Process sequence metadata after loading.

        Returns
        -------
        None

        Notes
        -----
        Sequence metadata currently requires no extra processing and therefore
        defers to the abstract base implementation.
        """
        return super().process_metadata()
    
    def update_metadata(self,key:str,field:str,value:str) -> None:
        """
        Replace a metadata value for a given row and field.

        Parameters
        ----------
        key : str
            Row identifier to update.
        field : str
            Column name to update.
        value : str
            Replacement value.

        Returns
        -------
        None
        """
        index=self.headers[field.upper()] - 1
        for i in self.data:
            if i == key:
                self.data[key][index]=value
    
    def get_headers(self) -> list[str]:
        """
        Return table headers in display order.

        Returns
        -------
        list
            Header names, with the key column annotated as ``KEY:<name>``.
        """
        headers=[]
        for key in self.headers:
            if self.headers[key]==self.key_index:
                headers.append(f'KEY:{key}')
            else:
                headers.append(key)
        return headers

    def get_values(self, field:str) -> dict[str, str]:
        """
        Retrieve non-empty values for a selected field.

        Parameters
        ----------
        field : str
            Field name to extract.

        Returns
        -------
        dict
            Mapping from row key to non-empty field value.
        """
        index=self.headers[field.upper()] - 1
        query={}
        for key in self.data:
            if self.data[key][index] != '':
                query[key]=self.data[key][index]
        return query
    
    def get_keys_per_notnull(self,field:str) -> list[str]:
        """
        Return keys with a non-empty value in a selected field.

        Parameters
        ----------
        field : str
            Field name to inspect.

        Returns
        -------
        list
            Row keys whose selected field is not empty.
        """
        keys=[key for key in self.data if self.data[key][self.headers[field.upper()]-1] != '']
        return keys
    
    def get_keys_per_2fields(self,field1:str,field2:str,value1:str,value2:str='') -> list[str]:
        """
        Return keys matching a pair of field constraints.

        Parameters
        ----------
        field1 : str
            First field name to inspect.
        field2 : str
            Second field name to inspect.
        value1 : str
            Required value in ``field1``.
        value2 : str, default ''
            Required value in ``field2``.

        Returns
        -------
        list
            Row keys satisfying both constraints.
        """
        keys=[key for key in self.data if (self.data[key][self.headers[field1.upper()]-1] == value1 and
              self.data[key][self.headers[field2.upper()]-1] == value2)]
        return keys
    
    def export_metadata(self,filename:str) -> None:
        """
        Export sequence metadata to a CSV file.

        Parameters
        ----------
        filename : str
            Output CSV path.

        Returns
        -------
        None
        """
        with open(filename,'w') as table:
            for header in self.headers.keys():
                table.write(f'{header}{self.delim}')
            table.write('\n')
            for key in self.data:
                table.write(f'{key}{self.delim}')
                for i in range(len(self.data[key])):
                    table.write(f'{self.data[key][i]}{self.delim}')
                table.write('\n')
    def get_strains_from_organism(self) -> dict[str, list[str]]:
        """
        Parse strain descriptors from the organism name field.

        Returns
        -------
        dict
            Mapping from accession to parsed host, location, strain number, and
            year values.
        """
        index=self.headers['ORGANISM_NAME']-1
        report={}
        pattern=re.compile(r'\(([^/]+)/([^/]+)/([^/]+)/([^/]+)\)',re.IGNORECASE)
        pattern2=re.compile(r'\(([^/]+)/([^/]+)/([^/]+)/([^/]+)/([^/]+)\)',re.IGNORECASE)
        for key in self.data:
            match=pattern.search(self.data[key][index])
            if match:
                location=match.group(2)
                strain_number = match.group(3)
                year = match.group(4)
                report[key]=['Homo sapiens',location,strain_number,year]
        for key in self.data:
            match=pattern2.search(self.data[key][index])
            if match:
                location=match.group(3)
                strain_number = match.group(4)
                year = match.group(5)
                host=match.group(2)
                report[key]=[host,location,strain_number,year]
            
        return report


#############################################################################################################################################################

class ClusterMetadata(MetadataTable):
    def __init__(self, filename: str) -> None:
        """
        Initialize a cluster metadata table backed by a TSV file.

        Parameters
        ----------
        filename : str
            Path to the cluster metadata TSV file.
        """
        super().__init__(filename, delim="\t" , key_index=0)
        self.process_metadata()
    
    def validate_metadata(self ,field:str) -> dict[str, Any]:
        """
        Identify clusters with multiple values in a selected field.

        Parameters
        ----------
        field : str
            Field name to inspect.

        Returns
        -------
        dict
            Mapping from cluster identifier to non-singular field values.
        """
        index=self.headers[field.upper()] - 1
        multi_value={}
        if index > 0 and index < len(self.headers)-1:
            for key in self.data:
                if len(self.data[key][index]) >1:
                    multi_value[key]=self.data[key][index]

        return multi_value
    
    def process_metadata(self) -> None:
        """
        Convert serialized dictionary fields into Python dictionaries.

        Returns
        -------
        None
        """
        for key in self.data:
             for i in range(1,len(self.data[key])):
                self.data[key][i]=ast.literal_eval(self.data[key][i])
    
    def update_metadata(self) -> None:
        """
        Update cluster metadata values.

        Returns
        -------
        None

        Notes
        -----
        Cluster metadata updates are not implemented and therefore defer to the
        abstract base implementation.
        """
        return super().update_metadata()
    
    def get_headers(self) -> list[str]:
        """
        Return table headers in display order.

        Returns
        -------
        list
            Header names, with the key column annotated as ``KEY:<name>``.
        """
        headers=[]
        for key in self.headers:
            if self.headers[key]==self.key_index:
                headers.append(f'KEY:{key}')
            else:
                headers.append(key)
        return headers
    def get_field(self,field:str) -> dict[str, Any]:
        """
        Return values from a selected field.

        Parameters
        ----------
        field : str
            Field name to extract.

        Returns
        -------
        dict
            Mapping derived from the selected field values.
        """
        index=self.headers[field.upper()] - 1
        col={key:value for key in self.data for value in self.data[key][index]}
        return col

    def get_H(self) -> dict[str, set[str]]:
        """
        Extract H subtype assignments from HA-containing clusters.

        Returns
        -------
        dict
            Mapping from cluster identifier to recovered H subtype tokens.

        Notes
        -----
        Cluster identifiers with mixed genotypes are written to
        ``mixedHclusters.json``.
        """
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
        with open('mixedHclusters.json','w') as pkl:
            json.dump(mixed,pkl)
        return clust_gen4
    
    def get_N(self) -> dict[str, set[str]]:
        """
        Extract N subtype assignments from NA-containing clusters.

        Returns
        -------
        dict
            Mapping from cluster identifier to recovered N subtype tokens.

        Notes
        -----
        Cluster identifiers with mixed genotypes are written to
        ``mixedNclusters.json``.
        """
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
        with open('mixedNclusters.json','w') as pkl:
            json.dump(mixed,pkl)
        return clust_gen6
    def query_clust(self,field:str,query:list[str]) -> None:
        """
        Print selected field values for a list of clusters.

        Parameters
        ----------
        field : str
            Field name to display.
        query : list
            Cluster identifiers to print.

        Returns
        -------
        None
        """
        index=self.headers[field.upper()] - 1
        for item in query:
            print(item,self.data[item][index])
    def convert_to_prop(self,field: str) -> dict[str, dict[str, float]]:
        """
        Convert cluster count dictionaries to proportions.

        Parameters
        ----------
        field : str
            Field name containing count dictionaries.

        Returns
        -------
        dict
            Mapping from cluster identifier to proportional values.
        """
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
    
    def export_metadata(self) -> None:
        """
        Export cluster metadata to disk.

        Returns
        -------
        None

        Notes
        -----
        Export behavior is not implemented here and defers to the abstract base
        implementation.
        """
        return super().export_metadata()

###############################################################################################################################################################

class BlastReportTable(MetadataTable):
    def __init__(self, filename:str) -> None:
        """
        Initialize a BLAST report table backed by a TSV file.

        Parameters
        ----------
        filename : str
            Path to the BLAST report TSV file.
        """
        super().__init__(filename, delim='\t',key_index=0)
    
    def load_data(self) -> dict[str, list[str]]:
        """
        Load the BLAST report and keep only best matches per sample.

        Returns
        -------
        dict
            Mapping from sample identifier to the selected report fields.
        """
        samples = []
        Sample=namedtuple('Sample',['SAMPLE_ID','BB_ACCESS','PERC_IDENTITY','SEGMENT','GENOTYPE','HOST'])
        with open(self.filename) as f:
            lines=f.readlines()
            for line in lines[1:]:
                fields=line.strip().split('\t')
                if fields[1]!='Unassigned':
                    samples.append(Sample(SAMPLE_ID=fields[0],BB_ACCESS=fields[1],PERC_IDENTITY=fields[2]\
                                      ,SEGMENT=fields[3],GENOTYPE=fields[4],HOST=fields[5]))
                else:
                    samples.append(Sample(SAMPLE_ID=fields[0],BB_ACCESS=fields[1],PERC_IDENTITY='NA',SEGMENT='NA',GENOTYPE='NA',HOST='NA'))
        to_keep={i.SAMPLE_ID:[] for i in samples}
        for s in samples:
            if to_keep[s.SAMPLE_ID]==[] and s.PERC_IDENTITY!='NA':
                best=0
                best=float(s.PERC_IDENTITY)
                to_keep[s.SAMPLE_ID].append(s.BB_ACCESS)
                to_keep[s.SAMPLE_ID].append(s.PERC_IDENTITY)
                to_keep[s.SAMPLE_ID].append(s.SEGMENT)
                to_keep[s.SAMPLE_ID].append(s.GENOTYPE)
                to_keep[s.SAMPLE_ID].append(s.HOST)
            elif to_keep[s.SAMPLE_ID]==[] and s.PERC_IDENTITY=='NA':
                to_keep[s.SAMPLE_ID]=[s.BB_ACCESS,s.PERC_IDENTITY,s.SEGMENT,s.GENOTYPE,s.HOST]
            else:
                if float(s.PERC_IDENTITY)>best:
                    best=float(s.PERC_IDENTITY)
                    to_keep[s.SAMPLE_ID]=[s.BB_ACCESS,s.PERC_IDENTITY,s.SEGMENT,s.GENOTYPE,s.HOST]
        return to_keep

    def process_metadata(self) -> None:
        """
        Process BLAST report metadata after loading.

        Returns
        -------
        None

        Notes
        -----
        BLAST reports are already normalized during loading, so this method
        defers to the abstract base implementation.
        """
        return super().process_metadata()
    
    def validate_metadata(self) -> None:
        """
        Validate BLAST report metadata.

        Returns
        -------
        None

        Notes
        -----
        Validation is not implemented for this class and therefore defers to
        the abstract base implementation.
        """
        return super().validate_metadata()
    
    def update_metadata(self) -> None:
        """
        Update BLAST report metadata values.

        Returns
        -------
        None

        Notes
        -----
        Updates are not implemented for this class and therefore defer to the
        abstract base implementation.
        """
        return super().update_metadata()
    
    def get_segment(self) -> dict[str, str]:
        """
        Return best-hit segment assignments.

        Returns
        -------
        dict
            Mapping from query accession to segment assignment.
        """
        id_seg={}
        field='segment'
        index=self.headers[field.upper()] - 1
        try:    
            for key in self.data:
                id_seg[key]=self.data[key][index]
        except IndexError:
            id_seg[key]='NA'
        return id_seg
    def get_genotype(self) -> dict[str, str]:
        """
        Return best-hit genotype assignments.

        Returns
        -------
        dict
            Mapping from query accession to genotype assignment.
        """
        id_gen={}
        field='genotype'
        index=self.headers[field.upper()] - 1
        try:    
            for key in self.data:
                id_gen[key]=self.data[key][index]
            
        except IndexError:
            id_gen[key]='NA'
        return id_gen
    def get_host(self) -> dict[str, str]:
        """
        Return best-hit host assignments.

        Returns
        -------
        dict
            Mapping from query accession to host assignment.
        """
        id_host={}
        field='host'
        index=self.headers[field.upper()] - 1
        try:    
            for key in self.data:
                id_host[key]=self.data[key][index]
        except IndexError:
            id_host[key]='NA'
        return id_host
    def export_metadata(self) -> None:
        """
        Export BLAST report metadata to disk.

        Returns
        -------
        None

        Notes
        -----
        Export behavior is not implemented here and defers to the abstract base
        implementation.
        """
        return super().export_metadata()
    
    def extract_h(self) -> dict[str, set[str]]:
        """
        Extract H subtype tokens from genotype assignments.

        Returns
        -------
        dict
            Mapping from sample identifier to a set of recovered H subtype
            tokens.
        """
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
    def extract_n(self) -> dict[str, set[str]]:
        """
        Extract N subtype tokens from genotype assignments.

        Returns
        -------
        dict
            Mapping from sample identifier to a set of recovered N subtype
            tokens.
        """
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


#################################################################################################################################################################

class ClusterReportTable(MetadataTable):
    def __init__(self, filename:str) -> None:
        """
        Initialize a cluster report table backed by a TSV file.

        Parameters
        ----------
        filename : str
            Path to the cluster report TSV file.
        """
        super().__init__(filename, delim='\t',key_index=0)
        self.process_metadata()

    def process_metadata(self) -> None:
        """
        Convert serialized report fields into Python dictionaries.

        Returns
        -------
        None
        """
        for key in self.data:
             for i in range(3,len(self.data[key])):
                self.data[key][i]=ast.literal_eval(self.data[key][i])

    def validate_metadata(self) -> None:
        """
        Validate cluster report metadata.

        Returns
        -------
        None

        Notes
        -----
        Validation is not implemented for this class and therefore defers to
        the abstract base implementation.
        """
        return super().validate_metadata()
    
    def update_metadata(self) -> None:
        """
        Update cluster report metadata values.

        Returns
        -------
        None

        Notes
        -----
        Updates are not implemented for this class and therefore defer to the
        abstract base implementation.
        """
        return super().update_metadata()
    
    def get_segment(self) -> dict[str, str | list[str]]:
        """
        Return segment assignments for clustered reports.

        Returns
        -------
        dict
            Mapping from query accession to one or more segment assignments.
        """
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
    def get_genotype(self) -> dict[str, str | list[str]]:
        """
        Return genotype assignments for clustered reports.

        Returns
        -------
        dict
            Mapping from query accession to one or more genotype assignments.
        """
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
    def get_host(self) -> dict[str, dict[str, int | float]]:
        """
        Return host assignments for clustered reports.

        Returns
        -------
        dict
            Mapping from query accession to host assignment dictionaries.
        """
        id_host={}
        field='hosts'
        index=self.headers[field.upper()] - 1
        for key in self.data:
            if self.data[key][0]!='NA':
                id_host[key]=self.data[key][index]
        
        return id_host
    def print_report(self) -> None:
        """
        Print the cluster report to standard output.

        Returns
        -------
        None
        """
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

    def extract_h(self) -> dict[str, set[str]]:
        """
        Extract H subtype tokens from clustered genotype assignments.

        Returns
        -------
        dict
            Mapping from sample identifier to a set of recovered H subtype
            tokens.
        """
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
    def extract_n(self) -> dict[str, set[str]]:
        """
        Extract N subtype tokens from clustered genotype assignments.

        Returns
        -------
        dict
            Mapping from sample identifier to a set of recovered N subtype
            tokens.
        """
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
    
    def convert_to_prop(self,field: str) -> dict[str, dict[str, float]]:
        """
        Convert clustered count dictionaries to proportions.

        Parameters
        ----------
        field : str
            Field name containing count dictionaries.

        Returns
        -------
        dict
            Mapping from sample identifier to proportional values.
        """
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
                        dic[i]=(dic[i]/tot)
                        dic[i]=round(dic[i],3)
                    perc[key]=dic
            except IndexError:
                continue
            except AttributeError:
                continue
        return perc


    def export_metadata(self) -> None:
        """
        Export cluster report metadata to disk.

        Returns
        -------
        None

        Notes
        -----
        Export behavior is not implemented here and defers to the abstract base
        implementation.
        """
        return super().export_metadata()

##############################################################################################################################################

class BlastClustReportTable(ClusterReportTable):
    def __init__(self, filename:str) -> None:
        """
        Initialize a hybrid BLAST/cluster report table.

        Parameters
        ----------
        filename : str
            Path to the combined report TSV file.
        """
        super().__init__(filename)
        self.process_metadata()
    def load_data(self) -> dict[str, list[Any]]:
        """
        Load the combined report and keep only best matches per sample.

        Returns
        -------
        dict
            Mapping from sample identifier to the selected report fields.
        """
        samples = []
        Sample=namedtuple('Sample',['SAMPLE_ID','CLUSTER_REP','CLUSTER','PERC_IDENTITY','SEGMENT','GENOTYPE','HOST'])
        with open(self.filename) as f:
            lines=f.readlines()
            for line in lines[1:]:
                fields=line.strip().split('\t')
                if fields[1]!='Unassigned':
                    samples.append(Sample(SAMPLE_ID=fields[0],CLUSTER_REP=fields[1],CLUSTER=fields[2],PERC_IDENTITY=fields[3]\
                                      ,SEGMENT=fields[4].replace(' ',''),GENOTYPE=fields[5].replace(' ',''),HOST=fields[6]))
                else:
                    samples.append(Sample(SAMPLE_ID=fields[0],CLUSTER_REP='NA',CLUSTER='NA',PERC_IDENTITY='NA',SEGMENT='NA',GENOTYPE='NA',HOST='NA'))
        to_keep={i.SAMPLE_ID:[] for i in samples}
        for s in samples:
            if to_keep[s.SAMPLE_ID]==[] and s.PERC_IDENTITY!='NA':
                best=0
                best=float(s.PERC_IDENTITY)
                to_keep[s.SAMPLE_ID].append(s.CLUSTER_REP)
                to_keep[s.SAMPLE_ID].append(s.CLUSTER)
                to_keep[s.SAMPLE_ID].append(s.PERC_IDENTITY)
                to_keep[s.SAMPLE_ID].append(s.SEGMENT)
                to_keep[s.SAMPLE_ID].append(s.GENOTYPE)
                to_keep[s.SAMPLE_ID].append(s.HOST)
            elif to_keep[s.SAMPLE_ID]==[] and s.PERC_IDENTITY =='NA':
                to_keep[s.SAMPLE_ID]=[s.CLUSTER_REP,s.CLUSTER,s.PERC_IDENTITY,s.SEGMENT,s.GENOTYPE,s.HOST]
            else:
                if float(s.PERC_IDENTITY)>best:
                    best=float(s.PERC_IDENTITY)
                    to_keep[s.SAMPLE_ID]=[s.CLUSTER_REP,s.CLUSTER,s.PERC_IDENTITY,s.SEGMENT,s.GENOTYPE,s.HOST]
        return to_keep
    def process_metadata(self) -> None:
        """
        Convert serialized combined-report fields into Python dictionaries.

        Returns
        -------
        None
        """
        for key in self.data:
             for i in range(4,len(self.data[key])):
                try:
                    if type(self.data[key][i])==str:
                        self.data[key][i]=ast.literal_eval(self.data[key][i])
                except NameError:
                    self.data[key][i]=str(self.data[key][i])
    def export_metadata(self) -> None:
        """
        Export combined BLAST/cluster metadata to disk.

        Returns
        -------
        None
        """
        return super().export_metadata()
    def validate_metadata(self) -> None:
        """
        Validate combined BLAST/cluster metadata.

        Returns
        -------
        Any
            Result returned by the parent implementation.
        """
        return super().validate_metadata()
    def extract_h(self) -> dict[str, set[str]]:
        """
        Extract H subtype assignments from genotype values.

        Returns
        -------
        dict
            Mapping from sample identifiers to recovered H subtype tokens.
        """
        return super().extract_h()
    def extract_n(self) -> dict[str, set[str]]:
        """
        Extract N subtype assignments from genotype values.

        Returns
        -------
        dict
            Mapping from sample identifiers to recovered N subtype tokens.
        """
        return super().extract_n()
    def print_report(self) -> None:
        """
        Print the combined report to standard output.

        Returns
        -------
        None
        """
        return super().print_report()
    def update_metadata(self) -> None:
        """
        Update combined BLAST/cluster metadata values.

        Returns
        -------
        None
        """
        return super().update_metadata()
    def convert_to_prop(self, field: str) -> dict[str, dict[str, float]]:
        """
        Convert count dictionaries in a report field to proportions.

        Parameters
        ----------
        field : str
            Field name containing count dictionaries.

        Returns
        -------
        dict
            Mapping from sample identifiers to proportional values.
        """
        return super().convert_to_prop(field)

##############################################################################################################################################
