#!/usr/bin/python3
import re
import subprocess
from Bio import SeqIO
import json
import random
import os
import shutil
from typing import Any, Iterable, Optional
from structures import segment_syns, iupac_nucleotides
import pandas as pd

def strains_get(filename: str) -> list[str]:
    """
    Collect strain descriptors from a multi-FASTA file.

    Parameters
    ----------
    filename : str
        Path to the FASTA file.

    Returns
    -------
    list
        Unique strain descriptors matched from FASTA headers.
    """
    strains=set()
    with open(filename,'r') as fasta:
        for line in fasta.readlines():
            if '>' in line:
                pattern = re.compile(r"\(.*\(.*\)\)", re.IGNORECASE)
                ex=re.search(pattern,line) 
                if ex is not None:
                    strains.add(ex.group())
    return list(strains)

def acess_get(strainlist: Iterable[str], filename: str) -> dict[str, list[str]]:
    """
    Collect accession numbers for selected strain names.

    Parameters
    ----------
    strainlist : list
        Strain names to search for.
    filename : str
        Path to the FASTA file.

    Returns
    -------
    dict
        Mapping from strain name to matching accession prefixes.
    """
    strain_seg={}
    for strain in strainlist:
        strain_seg[strain]=[]       #primes dict to recieve accessions
    with open(filename,'r') as fasta:
        for line in fasta.readlines():
            for key in strain_seg:
                if key in line:
                    strain_seg[key].append(line[1:12])

    return strain_seg

def seq_get(filename: str) -> dict[str, str]:
    """
    Parse a FASTA file into a header-to-sequence mapping.

    Parameters
    ----------
    filename : str
        Path to the FASTA file.

    Returns
    -------
    dict
        Dictionary mapping normalized FASTA headers to uppercase nucleotide
        sequences.

    Notes
    -----
    Sequences are normalized by removing ``-``, replacing ``X`` with ``N``,
    converting ``U`` to ``T``, and replacing non-IUPAC characters with ``N``.
    """
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



def seq_filter_get(filename: str,values: Iterable[str]) -> dict[str, str]:
    """
    Parse a FASTA file and keep only headers matching selected values.

    Parameters
    ----------
    filename : str
        Path to the FASTA file.
    values : iterable
        Substrings to match against FASTA headers.

    Returns
    -------
    dict
        Dictionary mapping matching FASTA headers to their sequences.
    """
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

def remove_seq_from_fasta(filename: str, output_path: str, access: list[str]) -> None:
    """
    Remove sequences from a FASTA file by accession.

    Parameters
    ----------
    filename : str
        Path to the input FASTA file.
    output_path : str
        Path to the output FASTA file.
    access : list
        Accession numbers to remove.

    Returns
    -------
    None
    """
    sequences = SeqIO.parse(filename,"fasta")
    filtered_sequences= (seq for seq in sequences if seq.id not in access)
    SeqIO.write(filtered_sequences, output_path, "fasta")

def sp_blastn(query: str,db: str,outputname: str = 'results',outputformat: str = '6',remote: bool = False,createdb: bool = True, maxtargetseqs: int = 1,silent: bool = False, numthreads: int = 1) -> None:
    """
    Run `blastn` against a local or remote nucleotide database.

    Parameters
    ----------
    query : str
        Path to the query FASTA file.
    db : str
        FASTA file used to build a local database, or an existing remote
        database name.
    outputname : str, default 'results'
        Output file stem for the BLAST results.
    outputformat : str, default '6'
        BLAST output format code.
    remote : bool, default False
        Whether to query a remote BLAST database.
    createdb : bool, default True
        Whether to create a local BLAST database from ``db``.
    maxtargetseqs : int, default 1
        Maximum number of aligned target sequences.
    silent : bool, default False
        Whether to suppress printing the result file contents.
    numthreads : int, default 1
        Number of CPU threads to pass to BLAST.

    Returns
    -------
    None
    """
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

def query_seqs(seqdict: dict[str, str],access: str) -> None:
    """
    Print FASTA headers containing a selected accession token.

    Parameters
    ----------
    seqdict : dict
        Mapping from FASTA headers to sequences.
    access : str
        Accession token to search for.

    Returns
    -------
    None
    """
    for key in seqdict:
        if access in key:
            print(key)

def to_subprocess(command: str) -> list[str]:
    """
    Split a shell-like command string into subprocess arguments.

    Parameters
    ----------
    command : str
        Command string to split on spaces.

    Returns
    -------
    list
        Tokenized subprocess argument list.
    """
    subprocess=command.split(' ')
    return subprocess

def sp_grep(command: str) -> None:
    """
    Run `grep` via `subprocess`.

    Parameters
    ----------
    command : str
        Grep arguments, with or without the leading ``grep`` token.

    Returns
    -------
    None
    """
    cmd=to_subprocess(command)
    greplist=['grep']
    if 'grep' not in cmd:
        greplist.append(cmd)
        subprocess.run(greplist)
    else:
        subprocess.run(cmd)

def filter_seq(seq_dict: dict[str, str], query: Iterable[str]) -> dict[str, str]:
    """
    Filter sequence dictionary entries by header substring.

    Parameters
    ----------
    seq_dict : dict
        Mapping from FASTA headers to sequences.
    query : iterable
        Substrings to match against each header.

    Returns
    -------
    dict
        Filtered dictionary containing matching entries.
    """
    filtered={}
    for key in seq_dict:
        for value in query:
            if key.upper().find(value.upper())!=-1:
                filtered[key]=seq_dict[key]
    return filtered

def concat_fasta(flist: list[str],outname: str) -> None:
    """
    Concatenate several FASTA files into one multi-FASTA file.

    Parameters
    ----------
    flist : list
        FASTA file paths to merge.
    outname : str
        Output file stem.

    Returns
    -------
    None
    """
    records=[]
    for i in flist:
        with open(i,'r') as file:
            fasta=file.readlines()
            name=''
            seq_lines=[]
            for line in fasta:
                if line.startswith('>'):
                    if name:
                        records.append((name, ''.join(seq_lines).upper()))
                    name=line.strip()
                    seq_lines=[]
                else:
                    seq_lines.append(line.strip())
            if name:
                records.append((name, ''.join(seq_lines).upper()))
    with open(f'{outname}.fasta','w') as output:
        for key, value in records:
            output.write(f'{key}\n{value}\n')

def dict_to_fasta(seqs: dict[str, str],outname:str) -> None:
    """
    Converts a dictionary into a FASTA file.

    Parameters
    ----------
    seqs : dict
        Dictionary with FASTA headers as keys and sequences as values.
    outname : str
        Output path WITHOUT the .fasta extension.

    Returns
    -------
    None
        Writes the output to '{outname}.fasta'.
    """
    with open(f'{outname}.fasta','w') as output:
        for key, value in seqs.items():
            output.write(f'{key}\n{value}\n')

def fasta_splitter(inputfile: str,outfile: str,batch_size: int = 8) -> None:
    """
    Split a multi-FASTA file into smaller FASTA batches.

    Parameters
    ----------
    inputfile : str
        Path to the input FASTA file.
    outfile : str
        Output file stem.
    batch_size : int, default 8
        Number of sequences per output FASTA file.

    Returns
    -------
    None
    """
    
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

def parse_clstr(filename: str, access_only: bool = True) -> dict[str, list[str]]:
    """
    Parse a CD-HIT ``.clstr`` file.

    Parameters
    ----------
    filename : str
        Path to the ``.clstr`` file.
    access_only : bool, default True
        Whether to keep only accession identifiers instead of the raw cluster
        lines.

    Returns
    -------
    dict
        Mapping from cluster label to member accessions or raw entries.
    """
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

def metadata_dict(filepath: str,col_list: list[int],sep: str = ',') -> dict[str, list[str]]:
    """
    Build a lookup dictionary from a delimited metadata table.

    Parameters
    ----------
    filepath : str
        Path to the metadata file.
    col_list : list
        Column indexes to extract for each accession.
    sep : str, default ","
        Field delimiter used in the metadata file.

    Returns
    -------
    dict
        Dictionary keyed by accession with a list of selected column values.
    """

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

def json_save(object: Any,filename: str) -> None:
    """
    Save a Python object as JSON.

    Parameters
    ----------
    object : Any
        Python object to serialize.
    filename : str
        Output filename stem without the ``.json`` extension.

    Returns
    -------
    None
    """
    
    with open(f'{filename}.json','w') as save:
        json.dump(object,save)
    
def json_load(filename: str) -> Any:
    """
    Load a Python object from a JSON file.

    Parameters
    ----------
    filename : str
        Path to the JSON file.

    Returns
    -------
    Any
        Deserialized Python object.
    """
    
    with open(filename,'r') as bin:
        object=json.load(bin)
    
    return object

def merge_fasta(input_files: list[str], output_file: str) -> None:
    """
    Merge multiple FASTA files into one.

    Parameters
    ----------
    input_files : list
        Paths to the input FASTA files.
    output_file : str
        Path to the merged output FASTA file.

    Returns
    -------
    None
    """
    with open(output_file, "w") as outfile:
        for input_file in input_files:
            with open(input_file, "r") as infile:
                for record in SeqIO.parse(infile, "fasta"):
                    SeqIO.write(record, outfile, "fasta")

    print(f"Merged {len(input_files)} files into {output_file}")

def extract_random_sequences(input:str, output:str, seqs:int ) -> None:
    """
    Extract a random subset of sequences from a FASTA file.

    Parameters
    ----------
    input : str
        Path to the input FASTA file.
    output : str
        Path to the output FASTA file.
    seqs : int
        Number of sequences to sample.

    Returns
    -------
    None
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

def concat_tabular(file_list: list[str],output_file: str) -> None:
    """
    Concatenate tabular files with duplicate row removal.

    Parameters
    ----------
    file_list : list
        Paths to tabular files with the same header structure.
    output_file : str
        Path to the concatenated output file.

    Returns
    -------
    None
    """
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

def headers_from_mult_fas(fasta_list: list[str],only_name: bool = False,out_list: bool = False) -> dict[str, list[str]] | list[str]:
    """
    Extract headers from one or more FASTA files.

    Parameters
    ----------
    fasta_list : list
        FASTA file paths to inspect.
    only_name : bool, default False
        Whether to return only normalized accession names.
    out_list : bool, default False
        Whether to return a list instead of a dictionary.

    Returns
    -------
    dict or list
        Header information in the requested representation.
    """
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

def mine_headers(header_dict: dict[str, str]) -> dict[str, list[str]]:
    """
    Infer genotype and segment values from FASTA headers.

    Parameters
    ----------
    header_dict : dict
        Mapping from accession identifiers to FASTA header strings.

    Returns
    -------
    dict
        Mapping from accession identifier to ``[genotype, segment]``.
    """
    
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

def input_missing_gen_seq(dataframe: pd.DataFrame, dictionary: dict,keydfname: str,value1dfname: str,value2: bool = False,value2dfname: Optional[str] = None) -> None:
    """
    Fill missing dataframe values from a lookup dictionary.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        Dataframe to update.
    dictionary : dict
        Lookup mapping keyed by dataframe identifiers.
    keydfname : str
        Name of the dataframe key column.
    value1dfname : str
        Name of the first dataframe column to fill.
    value2 : bool, default False
        Whether a second value should also be filled.
    value2dfname : str, optional
        Name of the second dataframe column to fill.

    Returns
    -------
    None
    """
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

def dict_to_csv(values_dict: dict[str, list[Any]],filename: str,header: Optional[list[str]] = [],sep: str = ',') -> None:
    """
    Write a dictionary of rows to a delimited text file.

    Parameters
    ----------
    values_dict : dict
        Mapping from row keys to row values.
    filename : str
        Output file path.
    header : list, default []
        Optional header row values.
    sep : str, default ','
        Output field delimiter.

    Returns
    -------
    None
    """
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

def mine_segment(dict: dict[str, str]) -> dict[str, str]:
    """
    Infer segment numbers from header text.

    Parameters
    ----------
    dict : dict
        Mapping from accession or sequence identifiers to header strings.

    Returns
    -------
    dict
        Mapping from identifier to inferred segment number.
    """
    mined={}
    for key,value in segment_syns.items():
        for access in dict.keys():
            for item in value:
                if item.upper() in dict[access].upper():
                    mined[access]=key
    return mined

def mine_genotype_H(list_of_genotypes: Iterable[str]) -> dict[str, list[str]]:
    """
    Extract HA subtype tokens from genotype labels.

    Parameters
    ----------
    list_of_genotypes : iterable
        Genotype strings such as ``H5N1`` or ``H5``.

    Returns
    -------
    dict
        Mapping from each input genotype string to the H subtypes recovered
        from it.
    """
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

def mine_single_HA(str: str) -> str:
    """
    Extract a single HA subtype token from a genotype string.

    Parameters
    ----------
    str : str
        Genotype string to inspect.

    Returns
    -------
    str
        Matching ``H`` subtype token, or ``"UD"`` when none is found.
    """
    part_h=r'H\d{1,2}'
    lookup_H=re.match(part_h,str)
    if lookup_H:
        return lookup_H.group()
    else:
        return 'UD'

def mine_single_NA(str: str) -> str:
    """
    Extract a single NA subtype token from a genotype string.

    Parameters
    ----------
    str : str
        Genotype string to inspect.

    Returns
    -------
    str
        Matching ``N`` subtype token, or ``"UD"`` when none is found.
    """
    part_n=r'N\d{1,2}'
    lookup_N=re.match(part_n,str)
    if lookup_N:
        return lookup_N.group()
    else:
        return 'UD'

def mine_genotype_N(list_of_genotypes: Iterable[str]) -> dict[str, list[str]]:
    """
    Extract NA subtype tokens from genotype labels.

    Parameters
    ----------
    list_of_genotypes : iterable
        Genotype strings such as ``H5N1`` or ``N1``.

    Returns
    -------
    dict
        Mapping from each input genotype string to the N subtypes recovered
        from it.
    """
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

def int_to_iupac(value: str | int) -> Optional[str]:
    """
    Convert a segment number into its IUPAC segment acronym.

    Parameters
    ----------
    value : str or int
        Segment number between 1 and 8.

    Returns
    -------
    str
        IUPAC acronym for the segment.
    """
    iupac={1:'PB2', 2: 'PB1', 3:'PA', 4:'HA', 5:'NP',6:'NA',7:'MP',8:'NS'}
    try:
        return iupac[int(value)]
    except KeyError:
        print(f"Invalid segment value: {value}, expected 1 to 8")

def convert_to_prop(dict: dict[str, float]) -> dict[str, float]:
    """
    Convert count values into proportions.

    Parameters
    ----------
    dict : dict
        Dictionary of count values.

    Returns
    -------
    dict
        Dictionary of proportional values.
    """
    prop={}
    total=sum(dict.values())
    for key in dict:
        prop[key]=dict[key]/total
    return prop


def get_ncbi_accessions(fasta_file: str) -> set[str]:
    """
    Extract NCBI accession numbers from FASTA headers.

    Parameters
    ----------
    fasta_file : str
        Path to the FASTA file.

    Returns
    -------
    set
        Extracted accession numbers.
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
