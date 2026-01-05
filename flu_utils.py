#!/usr/bin/python3
import re
import subprocess
from Bio import SeqIO
import json
import random
import os
import shutil
from structures import segment_syns, iupac_nucleotides
import pandas as pd

def strains_get(filename):
    '''Collects accession numbers from a multifasta
    Accepts: filename - a string for the path of the fasta file
    Returns: A list of accession numbers
    Raises: FileNotFoundError if the path is wrong'''
    strains=set()
    with open(filename,'r') as fasta:
        for line in fasta.readlines():
            if '>' in line:
                pattern = re.compile(r"\(.*\(.*\)\)", re.IGNORECASE)
                ex=re.search(pattern,line) 
                if ex is not None:
                    strains.add(ex.group())
    return list(strains)

def acess_get(strainlist, filename):
    '''Gets accession numbers from a list of strain names
    Accepts: strainlist - a list
    filename - a string (path to multifasta file)
    Returns: a dictionary with strain name as key and a list of accession
    numbers corresponding to that strain
    Raises: FileNotFoundError if the path is wrong'''
    strain_seg={}
    for strain in strainlist:
        strain_seg[strain]=[]       #primes dict to recieve accessions
    with open(filename,'r') as fasta:
        for line in fasta.readlines():
            for key in strain_seg:
                if key in line:
                    strain_seg[key].append(line[1:12])

    return strain_seg

def seq_get(filename):
    '''Parses fasta file and returns accession numbers and sequences
    Accepts: filename -  a string for the path of the fasta file
    Returns: a dictonary with accession number and sequence
    Raises:FileNotFoundError if file path is wrong'''
    seqs={}
    with open(filename,'r') as file:
        fasta=file.readlines()
    for i in range(len(fasta)):
        if '>' in fasta[i]:
            name=fasta[i].strip()
            name=name.replace(';','_')
            seqs[name]=''
        else:
            seqs[name]+=fasta[i].strip().upper().replace('-','').replace('X','N')
    for key in seqs:
        seqs[key]=seqs[key].replace('U','T')
        for nuc in seqs[key]:
            if nuc not in iupac_nucleotides:
                seqs[key]=seqs[key].replace(nuc,'N')
    return seqs



def seq_filter_get(filename,values):
    '''Parses fasta file and returns accession numbers and sequences
    Accepts: filename -  a string for the path of the fasta file
    values - an iterable with accession numbers or strains
    Returns: a dictonary with accession number and sequence
    Raises:FileNotFoundError if file path is wrong'''
    seqs={}
    seqsf={}
    with open(filename,'r') as file:
        fasta=file.readlines()
    for i in range(len(fasta)):
        if '>' in fasta[i]:
            name=fasta[i].strip()
            seqs[name]=''
        else:
            seqs[name]+=fasta[i].strip().upper()
    for key in seqs:
        for val in values:
            if val.upper() in key.upper():
                seqsf[key]=seqs[key]
    return seqsf

def remove_seq_from_fasta(filename, output_path, access: list):
    """
    Remove sequences from a FASTA file based on a list of accession numbers using BioPython.

    Parameters:
        filename (str): Path to the input FASTA file.
        output_path (str): Path to the output FASTA file.
        accession_to_remove (list): List of accession numbers to remove.

    Returns:
        None
    """
    sequences = SeqIO.parse(filename,"fasta")
    filtered_sequences= (seq for seq in sequences if seq.id not in access)
    SeqIO.write(filtered_sequences, output_path, "fasta")

def sp_blastn(query,db,outputname='results',outputformat='6',remote=False,createdb=True, maxtargetseqs=1,silent=False, numthreads=1):
    '''Creates a blast db from a nucleotide fasta file and runs a 
    blastn of a query sequence against the db
    Requires: ncbi-blast+ installed locally
    Accepts: query (string)- path to a .fasta file of the test sequence
    db (string)- filename of a multifasta file of sequences to target 
    that should be in the current working directory
    outputname (string)- name of the output file. (default- results)
    outputformat (string)- type of outformat (default- 6)
    remote(bool) - use with db=nt for blast against NCBI
    createdb(bool) -use to create a local blastDB from a fasta file
    maxtargetseqs(int) - default=1 max aligned seqs
    numthreads(int) - default=1 number of cpu threads
    Returns: .txt file and contents printed to stdout 
    '''
    dbname=db.split('.')[0]
    if createdb:
        subprocess.run(["makeblastdb", "-in", db, "-dbtype", "nucl", "-out",\
                     f"{dbname}"])
    if not remote:
        subprocess.run(["blastn", "-query", query, "-db", f"{dbname}", "-out",\
                     f"{outputname}.txt", "-outfmt", outputformat,"-subject_besthit","-max_target_seqs",f'{maxtargetseqs}',"-num_threads",f'{numthreads}'])
    else:
        subprocess.run(["blastn", "-query", query, "-db", f"{dbname}", "-out",\
                     f"{outputname}.txt", "-outfmt", outputformat,"-subject_besthit","-max_target_seqs",f'{maxtargetseqs}'])
    if not silent:
        with open(f"{outputname}.txt") as results_file:
            print(results_file.read())

def query_seqs(seqdict,access):
    '''Queries a sequence dictionary for the key with an accession number
    accepts: seqdict(dictionary) - dictionary of fasta header as key and 
    sequence as value
    access(str) - accession number
    prints header to stdout'''
    for key in seqdict:
        if access in key:
            print(key)

def to_subprocess(command):
    '''Transforms a string into a subprocess list
    Accepts: command - a string
    Returns: subprocess list - A list
    '''
    subprocess=command.split(' ')
    return subprocess

def sp_grep(command):
    '''Uses subprocess module to launch grep
    Accepts: command(str) command without grep
    Returns: prints results to stdout'''
    cmd=to_subprocess(command)
    greplist=['grep']
    if 'grep' not in cmd:
        greplist.append(cmd)
        subprocess.run(greplist)
    else:
        subprocess.run(cmd)

def filter_seq(seq_dict, query):
    '''Filters dictionary items by substring in key
    Accepts: seq_dict (dict)
    pattern(iterable)
    Returns: filtered dictionary
    '''
    filtered={}
    for key in seq_dict:
        for value in query:
            if key.upper().find(value.upper())!=-1:
                filtered[key]=seq_dict[key]
    return filtered

def concat_fasta(flist,outname):
    '''concatenates a list of fasta files and outputs a multifasta file
    Accepts: flist(list) - list of .fasta files
    outname(str) - name of the output fasta file
    Returns: fasta file'''
    seqs={}
    for i in flist:
        with open(i,'r') as file:
            fasta=file.readlines()
            for i in range(len(fasta)):
                if '>' in fasta[i]:
                    name=fasta[i].strip()
                    seqs[name]=''
                else:
                    seqs[name]+=fasta[i].strip().upper()
    with open(f'{outname}.fasta','w') as output:
        for key, value in seqs.items():
            output.write(f'{key}\n{value}\n')

def dict_to_fasta(seqs,outname:str):
    '''coverts a dictionary into a fasta file
    Accepts: seqs(dict) - dict with key as fasta header and values as seq
    outname(str) - name of the output fasta file
    Returns: fasta file'''
    with open(f'{outname}.fasta','w') as output:
        for key, value in seqs.items():
            output.write(f'{key}\n{value}\n')

def fasta_splitter(inputfile,outfile,batch_size=8):
    '''Requires Biopython
    Splits a multifasta file into smaller files with 
    batch_size sequences
    Accepts: inputfile(str) - path to fasta file
    outfile(str): name of output files 
    batch_size (int) number of sequences per fasta
    Returns: len(inputfile)//batch_size fasta files'''
    
    fasta = SeqIO.parse(inputfile, "fasta")
    counter=1
    seq_batch=[]

    for i, seq_record in enumerate(fasta, 1):
        seq_batch.append(seq_record)
        if i % batch_size == 0:
            out = f"{outfile}_{counter}.fasta"
            SeqIO.write(seq_batch, out, "fasta")
            print(f"Written {len(seq_batch)} sequences to {out}")
            seq_batch = []
            counter += 1
    if seq_batch:
        out = f"{outfile}_{counter}.fasta"
        SeqIO.write(seq_batch, out, "fasta")
        print(f"Written {len(seq_batch)} sequences to {out}")

def parse_clstr(filename, access_only=True):
    '''parses a .clstr file and retrieves the accessions within each cluster
    Accepts: filename(str) path to the clstr file
    access_only (bool) True for accession only, false for the entire row
    Returns: clusters(dict)'''
    clusters = {}
    current_cluster = None

    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('>Cluster'):
                current_cluster = line.strip()
                clusters[current_cluster] = []
            else:
                if not access_only:
                    clusters[current_cluster].append(line.strip())
                else:
                    line=line.split(' ')
                    line=line[1]
                    line=line.replace('...','')
                    line=line.replace('>','')
                    line=line.split('_|')[0]
                    clusters[current_cluster].append(line)
    return clusters

def metadata_dict(filepath,col_list,sep=','):
    '''generates a metadata dict from a metadata file with
    key=accession number; value=list of metadata columns
    cols, in flu_metadata.csv:
    0=Accession
    1=Organism_Name
    2=GenBank_RefSeq
    3=Assembly
    4=Release_Date
    5=Isolate
    6=Species
    7=Length
    8=Nuc_Completeness
    9=Genotype
    10=Segment
    11=Country
    12=Host
    13=Collection_Date
    
    Accepts: filepath(str) - path to metadata file
    col_list(list) - list of column indexes to mine
    Returns: dict - a metadata dictionary  '''

    #Primes dict
    dict={}

    #populates dict with metadata info defined by col_list
    with open(filepath,'r') as metadata:
        lines=metadata.readlines()
        for line in lines[1:]:
            line=line.strip()
            line=line.split(f'{sep}')
            dict[line[0]]=[]
            for col in col_list:
                dict[line[0]].append(line[col])

    return dict

def json_save(object,filename):
    """Saves an object into a json file
    Accepts: object(any) - Python object to store
    filename(str) - Path to the output filename w/o .json extension
    Outputs: .json file"""
    
    with open(f'{filename}.json','w') as save:
        json.dump(object,save)
    
def json_load(filename):
    """Loads an object from a json file
    Accepts: filename(str) - Path to the output filename (MUST INCLUDE .json)
    Outputs: python object"""
    
    with open(filename,'r') as bin:
        object=json.load(bin)
    
    return object

def merge_fasta(input_files, output_file):
    """
    Merges multiple FASTA files into one.
    REQUIRES: Biopython

    Accepts: input_files(list): paths to input FASTA files
    output_file (str): Path to the output FASTA file
    """
    with open(output_file, "w") as outfile:
        for input_file in input_files:
            with open(input_file, "r") as infile:
                for record in SeqIO.parse(infile, "fasta"):
                    SeqIO.write(record, outfile, "fasta")

    print(f"Merged {len(input_files)} files into {output_file}")

def extract_random_sequences(input:str, output:str, seqs:int ):
    """
    Extracts a subset of random sequences from a FASTA file.

    Accepts:input(str) - Path to the input FASTA file
    output(str) - Path to the output FASTA file
    seqs (int): Number of random sequences to extract
    """
    # Creates a list for all the seqs in the fasta file
    sequences = list(SeqIO.parse(input, "fasta"))

    # Check if the requested number of sequences is greater than available sequences
    if seqs > len(sequences):
        print(f"Only {len(sequences)} available.")
        seqs = len(sequences)

    # Select random sequences
    random_sequences = random.sample(sequences, seqs)
    # Write selected sequences to the output file

    with open(output, "w") as outfile:
        SeqIO.write(random_sequences, outfile, "fasta")
    
    print(f"Extracted {seqs} random sequences to {output}")

def concat_tabular(file_list,output_file):
    '''concatenates a list of tabular files with the same structure
    into a single tabular file, rejecting repeated lines
    Accepts: file_list(list) - list of paths to .txt tabular files
    output_file (str) - path to output file
    Outputs: concatentated tabular file (output_file_path_concat.txt)'''
    concatenated=set()
    #mining header
    header=open(file_list[0],'r').readlines()
    header=header[0]
    outlist=[]
    outlist.append(header)

    #opening files and adding the lines to the set
    for file in file_list:
        with open(file,'r') as tabular:
            lines=tabular.readlines()
            for line in lines:
                concatenated.add(line)

    #converting set to list and append values to outlist
    concatenated=list(concatenated)
    for line in concatenated:
        if line != header:
            outlist.append(line)
    
    #oututing the final report
    with open(output_file,'w') as output:
        for line in outlist:
            output.write(line)

def headers_from_mult_fas(fasta_list,only_name=False,out_list=False):
    '''queries fasta for headers and returns dict: 
    {Accesion:Name,...} if only_name==False
    {Accession:Name} if only_name == True
    Accepts: fasta(str - path to fasta file)
    only_name (bool): default=False
    Returns: headers (dict) '''
    if not only_name and not out_list:
        headers={}
    elif not only_name and out_list:
        headers=[]
    else:
        headers=[]
        
    for fasta in fasta_list:
        with open(fasta,'r') as file:
            lines=file.readlines()
            for line in lines:
                if not only_name and not out_list:   
                    if '>' in line:
                        header=line.strip().split('|')
                        header[0]=header[0].replace('>','')
                        header[0]=header[0].replace(' ','')
                        header[0]=header[0].replace(';','_')
                        headers[header[0]]=[]
                        for i in range(1,len(header)):
                            headers[header[0]].append(str(i))
                elif out_list and not only_name:
                    if '>' in line:
                        header=line.strip()
                        header=header.replace('>','')
                        header=header.replace(' ','_')
                        headers.append(header)
                else:
                    if '>' in line:
                        header=line.strip().split('|')
                        header[0]=header[0].replace('>','')
                        header[0]=header[0].replace(' ','')
                        header[0]=header[0].replace(';','_')
                        headers.append(header[0])
    return headers

def mine_headers(header_dict):
    '''Mines a dictionary from fasta headers to extract segment
    and genotype.
    Accepts: header_dict (dict) - a dictionary of fasta accession number: header.
    Returns: gen_seg (dict) - a dictionary of accession number: [genotype, segment].'''
    
    gen_seg = {}
    
    # Genotype regex patterns
    full_gen_patt = re.compile(r"H\d{1,2}N\d{1,2}", re.IGNORECASE)  # Full HxNx pattern
    gen_patt_H = re.compile(r"H\d{1,2}", re.IGNORECASE)  # Hx only pattern
    gen_patt_N = re.compile(r"N\d{1,2}", re.IGNORECASE)  # Nx only pattern

    for key in header_dict:
        # Initialize to 'UD' (Undefined)
        gen_seg[key] = ['UD', 'UD']
        
        # Retrieving genotype using regex
        ex = re.search(full_gen_patt, header_dict[key])
        
        if ex is not None:
            gen_seg[key][0] = ex.group()  # Full HxNx found
        else:
            # Try finding H or N individually if full HxNx not found
            ex1 = re.search(gen_patt_H, header_dict[key])
            ex2 = re.search(gen_patt_N, header_dict[key])
            
            if ex1 is not None and ex2 is not None:
                # Only set if both H and N are found individually
                gen_seg[key][0] = ex1.group() + ex2.group()
            elif ex1 is not None:
                gen_seg[key][0] = ex1.group()  # Only H found
            elif ex2 is not None:
                gen_seg[key][0] = ex2.group()  # Only N found
    
    # Segment matching
    for key in header_dict:
        for value, syn_list in segment_syns.items():
            for item in syn_list:
                gen_patt3 = re.compile(rf"\b{item}\b", re.IGNORECASE)  # Match whole words
                ex = re.search(gen_patt3, header_dict[key])
                
                if ex is not None:
                    gen_seg[key][1] = value
                    break  # Stop searching once a match is found
            if gen_seg[key][1] != 'UD':  # If a segment has been matched, stop outer loop too
                break
                
    return gen_seg

def input_missing_gen_seq(dataframe, dictionary,keydfname,value1dfname,value2=False,value2dfname=None):
    '''Inputs up to 2 missing values into a data frame (REQUIRES PANDAS)
    Accepts dictionary with key value that exist in the data_frame
    dataframe - a pandas dataframe
    keydfname(str) - name of the DF key column
    value1dfname(str) - name of the DF value 1 column
    value2(bool) - If there is a second value (default=False)
    value2dfname(str) - name of the DF value 2 column
    outputs updated dataframe'''
    for idx, row in dataframe.iterrows():
        key = row['Accession']
    
    # If the key exists in the dictionary
        if key in dictionary:
        # Check and update Value1 if it's NaN
            if pd.isna(row['Genotype']):
                dataframe.loc[idx, 'Genotype'] = dictionary[key][0]
            if value2:
                if pd.isna(row['Segment']):
                    dataframe.loc[idx, 'Segment'] = dictionary[key][1]

def dict_to_csv(values_dict,filename,header=[],sep=','):
    '''writes to a csv file a dictionary with the 
    row name as the row key and the other values as other 
    columns
    Acceptts: values_dict(a dictionary)
    headers (list) - list of str with header names 
    filename (str) - path and name of output file
    Outputs: csv file'''
    with open(f'{filename}','w') as csv:
        if header != []:    
            for i in header:
                csv.write(f'{i}{sep}')
            csv.write('\n')
        for key in values_dict:
            csv.write(f'{key}{sep}')
            for value in values_dict[key]:
                csv.write(f'{value}{sep}')
            csv.write('\n')

def mine_segment(dict):
    mined={}
    for key,value in segment_syns.items():
        for access in dict.keys():
            for item in value:
                if item.upper() in dict[access].upper():
                    mined[access]=key
def mine_genotype_H(list_of_genotypes):
    gen_patt=r'(H\d{1,2})(N\d{1,2})'
    part_h=r'H\d{1,2}'
    part_genos={gen:set() for gen in list_of_genotypes}
    for geno in list_of_genotypes:
        match= re.match(gen_patt,geno)
        if match:
            part_genos[geno].add(match.group(1)) #HA
        else:
            lookup_H=re.match(part_h,geno)
            if lookup_H:
                part_genos[geno].add(lookup_H.group())
    for key in part_genos:
        part_genos[key]=list(part_genos[key])
    return part_genos
def mine_single_HA(str):
    part_h=r'H\d{1,2}'
    lookup_H=re.match(part_h,str)
    if lookup_H:
        return lookup_H.group()
    else:
        return 'UD'
def mine_single_NA(str):
    part_n=r'N\d{1,2}'
    lookup_N=re.match(part_n,str)
    if lookup_N:
        return lookup_N.group()
    else:
        return 'UD'
def mine_genotype_N(list_of_genotypes):
    gen_patt=r'(H\d{1,2})(N\d{1,2})'
    part_n=r'N\d{1,2}'
    part_genos={gen:set() for gen in list_of_genotypes}
    for geno in list_of_genotypes:
        match= re.match(gen_patt,geno)
        if match:
            part_genos[geno].add(match.group(2)) #NA
        else:
            lookup_N=re.match(part_n,geno)
            if lookup_N:
                part_genos[geno].add(lookup_N.group())
    for key in part_genos:
        part_genos[key]=list(part_genos[key])
    return part_genos

def int_to_iupac(value):
    '''
    Converts a segment integer value into the IUPAC acronym
    Args: value(str or int) between 1 and 8
    Returns: str - IUPAC acronym for the segment 
    '''
    iupac={1:'PB2', 2: 'PB1', 3:'PA', 4:'HA', 5:'NP',6:'NA',7:'MP',8:'NS'}
    try:
        return iupac[int(value)]
    except KeyError:
        print(f"Invalid segment value: {value}, expected 1 to 8")

def convert_to_prop(dict):
    '''Converts a dictionary of counts into a dictionary of proportions
    Accepts: dict(dict) - dictionary of counts
    Returns: prop(dict) - dictionary of proportions'''
    prop={}
    total=sum(dict.values())
    for key in dict:
        prop[key]=dict[key]/total
    return prop


def get_ncbi_accessions(fasta_file):
    """
    Parses a multi-FASTA file and extracts NCBI accession numbers from headers.

    Args:
        fasta_file (str): Path to the FASTA file.

    Returns:
        set: A set of extracted NCBI accession numbers.
    """
    accessions = set()
    
    # NCBI accession number patterns
    accession_patterns = [
        r">([A-Z]{1,2}_\d+\.\d+)",  # RefSeq (e.g., "NM_001256789", "XP_123456")
        r">([A-Z]+\d+\.\d+)",  # GenBank (e.g., "ABC123456.1", "XYZ9876543.2")
        r">(\d+\.\d+)",        # Numeric-only (e.g., "123456.1")
    ]
    
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):  # FASTA headers start with '>'
                for pattern in accession_patterns:
                    match = re.search(pattern, line)
                    if match:
                        accessions.add(match.group(1))
                        break  # Stop after the first valid match
    
    return accessions




""" seg_dict=pkl_load('seg_dict.pkl')
mined_gen=mine_genotype(seg_dict)
mined_seg=mine_segment(seg_dict)"""
""" metadata_path='metadata/flu_metadata_(3.1)_.csv'
metadata=metadata_dict(metadata_path, [1,2,3,4,5,6,7,8,9,10,11,12,13])
mined_seg=pkl_load('seg_8.pkl')
for key in metadata:
    if metadata[key][9] not in ['1','2','3','4','5','6','7','8']:
        try:
            metadata[key][9]=mined_seg[key]
        except KeyError:
            continue
     """

""" outpath='metadata/flu_metadata_(3.1)_.csv'
header=['Accession','Organism_Name','GenBank_RefSeq','Assembly','Release_Date',
        'Isolate','Species','Length','Nuc_Completeness','Genotype','Segment',
        'Country','Host','Collection_Date']
dict_to_csv(metadata,header,outpath,sep=';') """