#!/home/joao/miniconda3/envs/fluid/bin/python3.10
#### IMPORT DEPENDENCIES
import os
import sys
import argparse
import configparser
import subprocess
import json
import glob
from flu_utils import seq_get,dict_to_fasta,headers_from_mult_fas,parse_clstr,pkl_load,\
    seq_filter_get,mine_genotype_H,mine_genotype_N,convert_to_prop, mine_single_HA,mine_single_NA
from metadata_utils import ClusterReportTable, ClusterMetadata, SequenceMetadata
from structures import flagdict
from copy import deepcopy
import pandas as pd
from gb_utils import fetch_genbank_list, append_genbank_from_list, get_gb, create_lookup,\
    write_gb,gb_to_fasta

### STEP FUNCTIONS
#### FASTA PREPROCESS
def seq_enum(counter: int) -> str:
    """
    Generate a sequence identifier string in the format 'Seq_{counter}'.

    Parameters:
    counter (int): The counter value to be used in the sequence identifier.

    Returns:
    str: A sequence identifier string in the format 'Seq_{counter}'.
    """
    return f'>Seq_{counter}'

def fasta_preprocess(filename: str,samples_path: str, min_seq_len: int, max_seq_len: int, flags:dict, verbose: bool=False) -> None:
    """
    Preprocess a FASTA file by removing sequences outside of a length range and saving the processed sequences and header mappings.

    Parameters:
    filename (str): The path to the FASTA file to be preprocessed.
    min_seq_len (int): The minimum sequence length threshold.
    max_seq_len (int): The maximum sequence length threshold.
    verbose (bool): Whether to print verbose output during preprocessing.

    Returns:
    Mappings (dict): A dictionary mapping sequence identifiers to sequence headers.
    """
    fasta_dict=seq_get(os.path.join(samples_path,filename))
    
    #removing length outliers
    
    to_pop=[]
    for key in fasta_dict:
        if len(fasta_dict[key]) < min_seq_len or len(fasta_dict[key]) > max_seq_len:
            to_pop.append(key)
    for item in to_pop:
        if "Fasta Preprocess" in flags:
            flags["Fasta Preprocess"]["Rejected Sequences"].append(item)
        fasta_dict.pop(item)
    if verbose:
        if len(to_pop)==1:
            print(f'{len(to_pop)} sequence was removed for being outside the length thresholds:')
        elif to_pop==[]:
            print('All sequences were within length thresholds')
        else:
            print(f'{len(to_pop)} sequences were removed for being outside the length thresholds:')
    
    #creating a dictionary for header mappings
    
    header_mapping={}
    counter=1
    for header in fasta_dict:
        key=seq_enum(counter)
        header_mapping[key]=header
        counter+=1
    
    #converting to processed dictionary
    
    mappings={}
    for key in header_mapping:
        mappings[key]=fasta_dict[header_mapping[key]]
    outname=filename.split('/')[-1]
    outname=outname.replace('.fasta','')
    outname_fasta='format_'+outname
    dict_to_fasta(mappings, os.path.join(samples_path,outname_fasta))
    return header_mapping

#### CD-HIT
def cd_hit_est_2d(filename: str, cluster_reps: str, output:str, identity: float, logs_p:str, verbose: bool=False ) -> None:
    '''
    Run CD-HIT-EST-2D to cluster sequences based on a given identity threshold.
    Parameters:
    filename (str): The path to the FASTA file to be clustered.
    cluster_reps (str): The path to the cluster representatives file.
    output (str): The path to the output file.
    identity (float): The sequence identity threshold for clustering.

    '''
    if not verbose:
        os.system(f'cd-hit-est-2d -i {cluster_reps} -i2 {filename} -o {output} -c {identity}\
                   >>{logs_p}{filename}_cd_hit.log 2>>{logs_p}{filename}_cd_hit.err')
    else:
        os.system(f'cd-hit-est-2d -i {cluster_reps} -i2 {filename} -o {output} -c {identity}')

#### CLUSTERING AND CLUSTER REPORT
def cluster_assign(report:str, headers_file:str,runsdir:str) -> dict:
    '''
    Assign sequences to clusters based on a CD-HIT-EST-2D output file and a FASTA headers file.
    Parameters:
    report (str): The path to the CD-HIT-EST-2D output file.
    headers_file (str): The path to the FASTA headers file.
    runsdir (str): The path to the runs directory.
    Returns:
    dict: A dictionary mapping cluster identifiers to assigned sequences.
    Outputs a .assign file with cluster assignments.
    '''
    #parsing clstr file
    clusters=parse_clstr(os.path.join(runsdir,report),access_only=False)
    report={}
    assigned=[]
    for cl in clusters:
        if clusters[cl][-1].split('\t')[0]!='0':
            report[cl]=[]
            for i in clusters[cl]:
                if "... *" not in i:
                    name=i.split(' ')[1]
                    name=name.replace('...','')
                    name=name.replace('>','')
                    name=name.replace('|','')
                    id=i.split(' ')[3]
                    id=id.replace('+/','')
                    id=id.replace('-/','')
                    assigned.append(name)
                    report[cl].append((name,id))
    #Extracting header information from the sample fasta
    headers=headers_from_mult_fas([headers_file],only_name=True)
    unassigned=[item.strip() for item in headers if item.strip() not in assigned]
    report['Unassigned']=unassigned
    sample_dict={}
    for header in headers:
        header=header.replace('|','')
        header=header.replace(' ','_')
        header=header.replace(';','_')
        header=header.replace('>','')  
        sample_dict[header]=''
    #Extracting identity and cluster information
    for key in report:   
        for i in range(len(report[key])):
            if type(report[key][i])==tuple:
                sample_dict[report[key][i][0]]=key
            else:
                sample_dict[report[key][i]]=f'>{key}'
    
    output_name=headers_file.split('/')[-1]
    output_name=output_name.replace('.fasta','.assign')
    outpath=os.path.join(runsdir,output_name)

    #Outputting to .assign file
    with open(f'{outpath}','w') as repo:
        repo.write('Cluster;Sample;%ID\n')
        for item in report:
            if item != 'Unassigned':
                for value in report[item]:
                    repo.write(f'{item};')
                    repo.write(f'{value[0]};')
                    repo.write(f'{value[1]}\n')
            else:
                for i in range(len(report['Unassigned'])):
                    repo.write(f'>Unassigned;{report["Unassigned"][i]}\n')
    return sample_dict

def cluster_compile(s_dict:str, cl_dict:dict, cl_assign:str, reports_path:str) -> None:
    '''
    Compile cluster assignments into a report file.
    Parameters:
    s_dict (str): The path to the sample dictionary file.
    cl_assign (str): The path to the cluster assignments file.
    cl_dict (dict): The cluster metadata dictionary.
    reports_path (str): The path to the reports directory.
    Outputs a compiled report file.
    '''
    #opening sample and cluster dictionaries and cluster assignments
    dict_sample=s_dict
    assign=open(f'{cl_assign}','r').readlines()[1:]
    dict_clust=cl_dict
    sample_assign={}
    
    #parsing cluster assignments
    for line in assign:
        line=line.strip()
        line=line.split(';')
        if len(line)>2:
            sample_assign[line[1]]=line[2]
        else:
            sample_assign[line[1]]='NA'
    compiler={}
    
    #adding to the information sample:cluster,rep,identity the seg,gen and host ranges
    for key in dict_sample:
        try:
            compiler[key]=[]
            compiler[key].append(sample_assign[key])
            compiler[key].append(dict_sample[key])
            if dict_sample[key]!='>Unassigned':
                for id in dict_clust[dict_sample[key]]:
                    compiler[key].append(id)
        except KeyError:
            print(f'Key: {key} not found')
            continue
    output_name=cl_assign.split('/')[-1]
    output_name=output_name.replace('.assign','')
    outpath=os.path.join(reports_path,f'{output_name}_clust_report.txt')
    
    #outputting report
    with open(f'{outpath}','w') as report:
        report.write('Sample_access\t%ID\tCluster\tCluster_rep\tGenotypes\tSegment\tHosts\n')
        for sample in compiler:
            if len(compiler[sample])>2:
                report.write(f'{sample}\t{compiler[sample][0]}\t\
                            {compiler[sample][1]}\t{compiler[sample][5]}\t\
                            {compiler[sample][2]}\t{compiler[sample][3]}\t{compiler[sample][4]}\n')
            elif len(compiler[sample])==0:
                continue  
            else:
                report.write(f'{sample}\t{compiler[sample][0]}\t{compiler[sample][1]}\n')
    print(f"Compiled report generated in {outpath}")

def cluster_miner(reports_p:str,filename:str,samples_p:str,flags:dict) -> dict:
    '''
    Mine cluster report for unassingned samples and 
    outputs them to BLAST
    Parameters:
    reports_p (str): The path to the reports directory.
    filename (str): The name of the sample file.
    samples_p (str): The path to the samples directory.
    runs_p (str): The path to the runs directory.
    Returns:
    dict: A dictionary mapping cluster identifiers to assigned sequences.
    Outputs a .fasta file for blast analysis.
    '''
    #mining the cluster report for unassigned samples
    clust_rep=ClusterReportTable(os.path.join(reports_p,f'format_{filename}_clust_report.txt'))
    clust_rep.convert_to_prop('genotypes')
    clust_rep.convert_to_prop('hosts')
    hs=clust_rep.extract_h()
    ns=clust_rep.extract_n()
    to_blast=[]
    for key in clust_rep.data:
        if clust_rep.data[key][0]=='NA': #completely unassigned
            to_blast.append(key)
        elif len(clust_rep.data[key][4]) > 1: #ambiguous in segment
            to_blast.append(key)
        elif clust_rep.data[key][4] == '4': #ambiguouis in subtype HA
            if len(hs[key]) > 1:
                to_blast.append(key)
        elif clust_rep.data[key][4] =='6': #ambiguouis in subtype NA
            if len(ns[key]) > 1:
                to_blast.append(key)
    to_report={}

    #creating a dictionary for the report
    for key in clust_rep.data:
        if key not in to_blast:
            to_report[key]=[]
            to_report[key].append(clust_rep.data[key][2].replace(' ',''))
            to_report[key].append(
                clust_rep.data[key][1].replace('                        >','>'))
            to_report[key].append(clust_rep.data[key][0].replace('%',''))
            to_report[key].append(clust_rep.data[key][4])
            to_report[key].append(clust_rep.data[key][3])
            to_report[key].append(clust_rep.data[key][5])
            to_report[key].append('CD-HIT')
    
    #swtting sequences to blast and updating flagsdict
    if to_blast != []:
        for i in to_blast:
            flags['CD-HIT']['Unclustered Sequences'].append(i)
        to_blast_dict=seq_filter_get(os.path.join(samples_p,f'format_{filename}.fasta'),to_blast)
        dict_to_fasta(to_blast_dict,os.path.join(samples_p,f'{filename}_to_blast'))
        flags['Master']['C_BLAST']=True
    else:
        flags['Master']['C_BLAST']=False
    return to_report

#### BLAST
def best_blast(runs_p:str,blast_report:str,thres_id:float=0.0) -> dict:
    '''
    Selects the best matches for blast reports and returns them as a dict
    Parameters:
    runs_p (str): path to the runs directory
    blast_report (str): name of the blast output file
    thres_id (float) - default = 0.0: minimum threshold of identity to accept a sequence
    Returns:
    best_blast(dict): dictionary with {sequence_name:[representative, identity]}}
    '''
    blast_df=pd.read_table(os.path.join(runs_p,blast_report),header=None,\
                           names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",\
                                   "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
    df_sorted = blast_df.sort_values(by=["qseqid", "evalue", "pident"], ascending=[True, True, False])
    best_hits = df_sorted.groupby("qseqid").first().reset_index()
    best_hits=best_hits[['qseqid','sseqid','pident']]
    best_hits=best_hits[best_hits['pident']>=thres_id]
    out_dict={}
    for row in best_hits.itertuples(index=False):
        out_dict[row.qseqid]=[row.sseqid,row.pident]
    return out_dict

def bclust(metadata_p:str,metadata_f:str,samples_p:str,blast_p:str,blast_db:str,runs_p:str,filename:str,flags:dict,fasta:str='to_blast.fasta') -> dict:
    '''
    Run BLAST analysis on unassigned samples.
    Parameters:
    metadata_p (str): The path to the metadata directory.
    samples_p (str): The path to the samples directory.
    blast_p (str): The path to the BLAST database directory.
    runs_p (str): The path to the runs directory.
    filename (str): The name of the sample file.
    fasta (str): The name of the FASTA file to be analyzed
    Outputs a .fasta file for BLAST analysis.
    '''
    access=headers_from_mult_fas([os.path.join(samples_p,fasta)],only_name=True)
    queries=[key for key in access]
    os.system(f'blastn -db {os.path.join(blast_p,blast_db)} -query {os.path.join(samples_p,fasta)}\
               -out {os.path.join(runs_p,f"{filename}_bcrun.txt")} -outfmt 6 -num_threads 2 -max_target_seqs 7')
    b_tab=best_blast(runs_p,f"{filename}_bcrun.txt",thres_id=90.0)
    pd.set_option('display.max_colwidth', None)
    metadata=pd.read_table(os.path.join(metadata_p,metadata_f),index_col=False)
    assigned=set(b_tab.keys())
    queries=set(queries)
    to_reblast=queries-assigned
    for key in b_tab:
        b_tab[key].append(metadata['Cluster'][metadata['Representative']==b_tab[key][0]].to_string(index=False,min_rows=None,max_rows=None))
        b_tab[key].append(metadata['Segments'][metadata['Representative']==b_tab[key][0]].to_string(index=False,min_rows=None,max_rows=None))
        b_tab[key][3]=eval(b_tab[key][3])
        b_tab[key].append(metadata['Genotypes'][metadata['Representative']==b_tab[key][0]].to_string(index=False,min_rows=None,max_rows=None))
        b_tab[key][4]=eval(b_tab[key][4])
        b_tab[key].append(metadata['Hosts'][metadata['Representative']==b_tab[key][0]].to_string(index=False,min_rows=None,max_rows=None))
        b_tab[key][5]=eval(b_tab[key][5])
    if list(to_reblast) != []:
        for i in to_reblast:
            flags['BLAST']['Sequences unassigned against cluster representatives'].append(i)
        to_reblast_dict=seq_filter_get(os.path.join(samples_p,fasta),to_reblast)
        dict_to_fasta(to_reblast_dict,os.path.join(samples_p,f'{filename.replace(".fasta","")}_to_reblast'))
        flags['Master']['L_BLAST']=True
    else:
        flags['Master']['L_BLAST']=False
    return b_tab

def reblast(metadata_p:str,metadata_f:str,samples_p:str,blast_p:str,blast_db:str,runs_p:str,filename:str,fasta:str='to_reblast.fasta') -> dict:

    '''
    Re-run BLAST analysis on unassigned samples.
    Parameters:
    metadata_p (str): The path to the metadata directory.
    samples_p (str): The path to the samples directory.
    blast_p (str): The path to the BLAST database directory.
    runs_p (str): The path to the runs directory.
    filename (str): The name of the sample file.
    fasta (str): The name of the FASTA file to be analyzed

    '''
    os.system(f'blastn -db {os.path.join(blast_p,blast_db)} -query {os.path.join(samples_p,fasta)}\
               -out {os.path.join(runs_p,f"{filename}_brun.txt")} -outfmt 6 -num_threads 2 -max_target_seqs 7')
    report=best_blast(runs_p,f"{filename}_brun.txt")
    metadata=pd.read_csv(os.path.join(metadata_p,metadata_f),index_col=False)
    for key in report:
        report[key].append(metadata['SEGMENT'][metadata['ACCESSION']==report[key][0]].to_string(index=False))
        report[key].append(metadata['GENOTYPE'][metadata['ACCESSION']==report[key][0]].to_string(index=False))
        report[key].append(metadata['HOST'][metadata['ACCESSION']==report[key][0]].to_string(index=False))
    return report

#### REPORT COMPILER
def report_compiler(clust_dict:dict,samples_p:str,filename:str,mappings_dict:dict,reports_p:str,flags:dict,bclust_dict:dict={},blast_dict:dict={}) -> None:
    '''
    Compile BLAST and cluster reports into a final report.
    Parameters:
    clust_dict (dict): The cluster report dictionary.
    samples_p (str): Path to samples file
    bclust_dict (dict): The BLAST cluster report dictionary.
    blast_dict (dict): The BLAST report dictionary.
    filename (str): The name of the sample file.
    mappings_dict (dict): The mappings dictionary.
    reports_p (str): The path to the reports directory.
    Outputs a final report file.
    '''
    seqs=headers_from_mult_fas([os.path.join(samples_p,f'format_{filename}')],only_name=True,out_list=True)
    to_report=clust_dict
    if bclust_dict:
        for key in bclust_dict:
            to_report[key]=[]
            to_report[key].append(bclust_dict[key][0])
            to_report[key].append(bclust_dict[key][2])
            to_report[key].append(bclust_dict[key][1])
            to_report[key].append(bclust_dict[key][3])
            to_report[key].append(convert_to_prop(bclust_dict[key][4]))
            to_report[key].append(convert_to_prop(bclust_dict[key][5]))
            to_report[key].append('C-BLAST')
    if blast_dict:
        for key in blast_dict:
            to_report[key]=[]
            to_report[key].append(blast_dict[key][0])
            to_report[key].append('UNCLUSTERED')
            to_report[key].append(blast_dict[key][1])
            to_report[key].append(blast_dict[key][2].replace(' ',''))
            to_report[key].append(blast_dict[key][3].replace(' ',''))
            to_report[key].append(blast_dict[key][4])
            to_report[key].append('L-BLAST')
    assigned=set(to_report.keys())
    seqs=set(seqs)
    to_remote=seqs-assigned
    mappings=mappings_dict
    if list(to_remote)!=[]:
        for seq in to_remote:
            to_report[seq]=['Unassigned','NA','NA','NA','NA','NA','NA']
    with open(f'{os.path.join(reports_p,filename.replace(".fasta",""))}_ID_Report.txt','w') as report:
        report.write('Sample_name\tRepresentative\tCluster\t%ID\tSegment\tGenotype\tHost\tAssigned_by\n'.upper())
        for key in to_report:
            mapped=mappings[f'>{key}']
            report.write(f'{mapped}\t{to_report[key][0]}\t\
                         {to_report[key][1]}\t{to_report[key][2]}\t{to_report[key][3]}\t\
                            {to_report[key][4]}\t{to_report[key][5]}\t{to_report[key][6]}\n')
    if len(list(to_remote))>0:
        flags['BLAST']['Sequences unassigned against local database'].extend(to_remote)
    print(f'Report generated in {os.path.join(reports_p,filename.replace(".fasta",""))}_ID_Report.txt')

#### REDIRECTOR
def redirector(report:str,flags:dict,filename:str,mappings:dict,samples_p:str,reports_p:str,force_flumut:bool,force_genin:bool,force_getref:bool,mode,single_sample=True) -> None:
    '''
    Redirects samples to aditional post-identification tools.
    Linchpin function that populates the 'Sample' and 'Final_report' dicts of flagsdict
    Parameters:
    report (str): The path to the final report file.
    flags (dict): The flags dictionary.

    '''
    #opening report and extracting data
    report=pd.read_table(os.path.join(reports_p,report), index_col=False)
    report=report.dropna()
    report['SEGMENT']=report['SEGMENT'].apply(lambda x: eval(x))
    report['SEGMENT']=report['SEGMENT'].apply(lambda x: list(x.keys())[0] if type(x)==dict else int(x))
    Segments=report['SEGMENT'].to_list()
    Genotypes=report['GENOTYPE'].to_list()
    Seq_names=report['SAMPLE_NAME'].to_list()
    #Type conversions and dictionary creations
    for i in range(len(Genotypes)):
        try:
            Genotypes[i]=eval(Genotypes[i])
        except NameError:
            continue
        Genotypes[i]=list(Genotypes[i].keys()) if\
            type(Genotypes[i])==dict else Genotypes[i]
    Segments=[int(i) for i in Segments]
    Seq_gen=dict(zip(Seq_names,Genotypes))
    Seq_seg=dict(zip(Seq_names,Segments))
    Seq_ref={k:v for k in report['SAMPLE_NAME'] for v in\
        report['REPRESENTATIVE'][report['SAMPLE_NAME']==k]}
    Seq_HA=''
    Seq_NA=''
    remap={v:k for k,v in mappings.items()}
    for key in Seq_seg:
        if Seq_seg[key]==4:
            Seq_HA=key
            Gen_HA=Seq_gen[Seq_HA]
            if type(Gen_HA)==list:
                for i in Gen_HA:
                    HA=mine_single_HA(i)
                    if HA != 'UD' and mode=='consensus':
                        if HA in flags['Final Report']['Sequences for NextClade']:
                            flags['Final Report']['Sequences for NextClade'][HA].append(remap[key])
            elif type(Gen_HA)==str:
                HA=mine_single_HA(Gen_HA)
                if HA != 'UD' and mode=='consensus':
                    if HA in flags['Final Report']['Sequences for NextClade']:
                        flags['Final Report']['Sequences for NextClade'][HA].append(remap[key])
            flags['Sample']['HA'].append(remap[key])
            flags['Sample']['HA_ref'].append(Seq_ref[key])
        elif Seq_seg[key]==6:
            Seq_NA=key
            Gen_NA=Seq_gen[Seq_NA]
            flags['Sample']['NA'].append(remap[key])
            flags['Sample']['NA_ref'].append(Seq_ref[key])
        elif Seq_seg[key]==1:
            flags['Sample']['PB2'].append(remap[key])
            flags['Sample']['PB2_ref'].append(Seq_ref[key])
        elif Seq_seg[key]==2:
            flags['Sample']['PB1'].append(remap[key])
            flags['Sample']['PB1_ref'].append(Seq_ref[key])
        elif Seq_seg[key]==3:
            flags['Sample']['PA'].append(remap[key])
            flags['Sample']['PA_ref'].append(Seq_ref[key])
        elif Seq_seg[key]==5:
            flags['Sample']['NP'].append(remap[key])
            flags['Sample']['NP_ref'].append(Seq_ref[key])
        elif Seq_seg[key]==7:
            flags['Sample']['MP'].append(remap[key])
            flags['Sample']['MP_ref'].append(Seq_ref[key])
        elif Seq_seg[key]==8:
            flags['Sample']['NS'].append(remap[key])
            flags['Sample']['NS_ref'].append(Seq_ref[key])
    if force_getref:
        flags['Final Report']['Get References'].extend(list(Seq_ref.values()))
    #Mining genotypes
    Hx=set()
    Nx=set()
    if single_sample:
        if 4 in Seq_seg.values():
            Gen_HA_mine=mine_genotype_H(Gen_HA)
            for key in Gen_HA_mine:
                try:
                    Hx.add(Gen_HA_mine[key][0])
                except IndexError:
                    continue
        if 6 in Seq_seg.values():
            Gen_NA_mine=mine_genotype_N(Gen_NA)
            for key in Gen_NA_mine:
                try:
                    if len(Gen_NA_mine[key])>1:
                        Nx.add(Gen_NA_mine[key][1])
                    else:
                        Nx.add(Gen_NA_mine[key][0])
                except IndexError:
                    continue
           
    if len(Hx)==1:
        Gen_HA=list(Hx)[0]
        flags['Sample']['H_gen']=Gen_HA
    if len(Nx)==1:
        Gen_NA=list(Nx)[0]
        flags['Sample']['N_gen']=Gen_NA
    if flags['Sample']['H_gen'] and flags['Sample']['N_gen']:
        flags["Sample"]["Genotype"]=f'{Gen_HA}{Gen_NA}'
    elif flags['Sample']['H_gen'] or flags['Sample']['N_gen']:
        flags["Sample"]["Genotype"]=f'{Gen_HA}' if flags['Sample']['H_gen'] else f'{Gen_NA}'
    else:
        flags["Sample"]["Genotype"]='Undetermined'
    #Directing samples to Final Reports
    if mode == 'consensus':
        if single_sample==False:
            force_flumut==True
            force_genin==True    
        if flags['Sample']['H_gen']=='H5' or force_flumut:
            for i in Seq_names:
                flags['Final Report']['Sequences for FluMut'].append(remap[i])
        if flags['Sample']['H_gen']=='H5' or force_genin:
            for i in Seq_names:
                flags['Final Report']['Sequences for GenIn'].append(remap[i])

    #Opening formatted fasta 
        if single_sample:
            formatted=seq_get(os.path.join(samples_p,f'format_{filename}'))
            flags['Sample']['HA_len']=len(formatted[flags['Sample']['HA'][0]])
            flags['Sample']['NA_len']=len(formatted[flags['Sample']['NA'][0]])
            flags['Sample']['PB2_len']=len(formatted[flags['Sample']['PB2'][0]])
            flags['Sample']['PB1_len']=len(formatted[flags['Sample']['PB1'][0]])
            flags['Sample']['PA_len']=len(formatted[flags['Sample']['PA'][0]])
            flags['Sample']['NP_len']=len(formatted[flags['Sample']['NP'][0]])
            flags['Sample']['MP_len']=len(formatted[flags['Sample']['MP'][0]])
            flags['Sample']['NS_len']=len(formatted[flags['Sample']['NS'][0]])
    elif mode=='contig':
        flags['Final Report']['Get References']=[value for value in Seq_ref.values()]
        if force_flumut:
            for i in Seq_names:
                flags['Final Report']['Sequences for FluMut'].append(remap[i])
        if force_genin:
            for i in Seq_names:
                flags['Final Report']['Sequences for GenIn'].append(remap[i])


#### FLUMUT PATH
def update_flumut_db():
    os.system('flumut --update')

def conform_to_flumut(flagdict: dict, sample_path: str, filename: str) -> None:
    """
    Adjusts sequence headers to conform to the flumut format and writes the adjusted sequences to a new FASTA file.
    Parameters:
    flagdict (dict): A dictionary containing flags and sequence information, including which sequences are designated for flumut processing.
    sample_path (str): The directory path where the input FASTA file is located.
    filename (str): The name of the FASTA file to be processed.

    Returns:
    None: This function does not return a value. It writes the adjusted sequences to a new FASTA file named 'to_flumut' in the specified sample path.
    """
    seq_seg = {}
    for key in flagdict['Sample']:
        if key in ('PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS'):
            for i in flagdict['Sample'][key]:
                if i in flagdict['Final Report']["Sequences for FluMut"]:
                    seq_seg[i] = key
    flt_fasta = seq_get(os.path.join(sample_path, f'format_{filename}'))
    outdict = {}
    for key in flt_fasta:
        if key in seq_seg:
            outdict[f'{key}|{seq_seg[key]}'] = flt_fasta[key]
    dict_to_fasta(outdict, os.path.join(sample_path, f'{filename.replace(".fasta","")}_to_flumut'))
        
def run_flumut(reports_p: str, samples_p: str, filename: str, regex="(.+)\|(.+)"):
    """
    Executes the flumut tool to analyze sequence data and generate mutation reports.

    Parameters:
    reports_p (str): The path to the directory where the report files will be stored.
    samples_p (str): The path to the directory containing the sample FASTA file.
    filename (str): The base name of the sample file, used to generate report filenames.
    regex (str): A regular expression pattern used by flumut for parsing sequence headers. Defaults to "(.+)\|(.+)".

    Returns:
    None: This function does not return a value. It executes a system command to run flumut and generate output files.
    """
    print("Running flumut")
    os.system(f'flumut -m {os.path.join(reports_p,f"{filename}_markers.tsv")} -M {os.path.join(reports_p,f"{filename}_mutations.tsv")} -l {os.path.join(reports_p,f"{filename}_literature.tsv")} {os.path.join(samples_p,f"{filename}_to_flumut.fasta")} -n "{regex}"')
def remap_flumut_report(reports_dir: str, mappings: dict, filename: str) -> None:
    """
    Remaps flumut report to the sample names.
    Parameters:
    reports_dir (str): The directory path where the reports are stored.
    mappings (dict): A dictionary containing the mapping of original sequence names to new names.
    filename (str): The name of the samples file

    Returns:
    None. The function modifies the flumut report files in-place by replacing the original sequence names with the new names.
    """
    # Create a mapping dictionary with the original sequence names replaced by the new names
    mapping_dict = {key.replace('>', ''): value for key, value in mappings.items()}

    # Read the markers report file and replace the original sequence names with the new names
    mk_report = pd.read_table(os.path.join(reports_dir,f'{filename}_markers.tsv'))
    mk_report = mk_report.replace(mapping_dict)

    # Read the mutations report file and replace the original sequence names with the new names
    mut_report = pd.read_table(os.path.join(reports_dir,f'{filename}_mutations.tsv'))
    mut_report = mut_report.replace(mapping_dict)

    # Save the modified markers report file
    mk_report.to_csv(os.path.join(reports_dir,f'{filename}_markers.tsv'), sep='\t', index=False)

    # Save the modified mutations report file
    mut_report.to_csv(os.path.join(reports_dir,f'{filename}_mutations.tsv'), sep='\t', index=False)

#### GET REFERENCE PATH
def get_reference(listref,references_p,db_p,filename,email=None, update_local_db=False):
    """
    Retrieves reference genomes from GenBank and stores them locally.
    Parameters:
    listref (list): A list of reference genome names.
    references_p (str): The path to the directory where the reference genomes will be stored."""
    references=set(listref)
    found=[]
    missing=[]

    #checking if local_db has entries
    records=create_lookup(os.path.join(references_p,db_p))
    for ref in references:
        if ref in records:
            found.append(ref)
        else:
            missing.append(ref)
    #getting missing records
    remote=fetch_genbank_list(missing,email)
    #getting local records
    local=get_gb(os.path.join(references_p,db_p), found)
    #creating output files
    output=remote+local
    #writing output files (.gb and .fasta)
    outname=f'{filename.replace(".fasta","")}_references'
    write_gb(output,os.path.join(references_p,f'{outname}.gb'))
    gb_to_fasta(os.path.join(references_p,f'{outname}.gb'),os.path.join(references_p,f'{outname}.fasta'))

    if update_local_db:
        append_genbank_from_list(os.path.join(references_p,db_p),remote)
#### NEXTCLADE PATH
def fasta_to_nextclade(flagsdict, samples_p, filename):
    to_nextclade={'H1':[],'H3':[],'H5':[]}
    for key in flagsdict['Final Report']['Sequences for NextClade']:
        to_nextclade[key].extend(flagsdict['Final Report']['Sequences for NextClade'][key])
        for key in to_nextclade:
            seqs=[]
            for access in to_nextclade[key]:
                seqs.append(access)
            flt_fasta = seq_get(os.path.join(samples_p, f'format_{filename}'))
            outdict = {}
            for sample in seqs:
                outdict[sample]=flt_fasta[sample]
            dict_to_fasta(outdict, os.path.join(samples_p, f'{filename.replace(".fasta","")}_to_nextclade_{key}'))
def run_nextclade(fasta,flagsdict,reports_p,samples_p,filename):
    dict_builds_broad={'H1':'flu_h1n1pdm_ha_broad','H3':'flu_h3n2_ha_broad','H5':'community/moncla-lab/iav-h5/ha/all-clades'}
    for key in flagsdict['Final Report']['Sequences for NextClade']:
        if flagsdict['Final Report']['Sequences for NextClade'][key] != []:
            fasta_gen=key    
            os.system(f'nextclade run -d {dict_builds_broad[fasta_gen]} --output-tsv {os.path.join(reports_p,f"{filename}_{fasta_gen}_nextclade.tsv")} {os.path.join(samples_p,f"{fasta}_{fasta_gen}.fasta")}')
def remap_nextclade(reports_p,nextclade_report,mappings_dict):
    """
    Remaps NextClade report to the sample names.
    Parameters:
    reports_p (str): The directory path where the reports are stored.
    nextclade_report (str): The path to the NextClade report file.
    mappings_dict (dict): A dictionary containing the mapping of original sequence names to new names.
    """
    nextclade_rep=pd.read_table(os.path.join(reports_p,nextclade_report),index_col=False)
    nextclade_rep['seqName']=nextclade_rep['seqName'].apply(lambda x: mappings_dict[f'>{x}'])
    nextclade_rep.to_csv(os.path.join(reports_p,nextclade_report), sep='\t', index=False)
    
def parser():
    parser=argparse.ArgumentParser(prog='AFluID',description='AFluID: Automated Influenza Identification Pipeline')
    parser.add_argument('-c','--config',type=str,help='Configuration file path',default='config.ini',required=False)
    parser.add_argument('-f','--filename',type=str,help='FASTA file name',required=True)
    parser.add_argument('-m','--mode',type=str,help='Mode of operation', choices=('contig','consensus'), required=True)
    parser.add_argument('-ff','--force',help='Force run of aditional tools',nargs='*',choices=('flumut','genin','getref'))
    parser.add_argument('-fdb','--update_flumut_db',type=str,help='turn off auto-update for flumut db', choices=('on','off'), default='on',required=False)
    parser.add_argument('-ss','--single_sample',type=str,help='Single sample mode',default='on',choices=('on','off'),required=False)
    parser.add_argument('-off','--turn_off',help='Turn Off additional analysis tools',nargs='*',choices=('flumut','genin','nextclade','getref'))
    parser.add_argument('-rm','--remove_previous',type=bool,help='Remove previous files',default=True,choices=(True,False),required=False)
    parser.add_argument('-Ml','--max_length',type=int,help='Maximum sequence length',required=False)
    parser.add_argument('-ml','--min_length',type=int,help='Minimum sequence length',required=False)
    args=parser.parse_args()
    return args

def main(flagdict=flagdict):
    ''' Main function for the pipeline'''
   
    ####LOAD FLAG DICT
    flags=deepcopy(flagdict)
    args=parser()
    filename=args.filename
    config_file=args.config
    mode=args.mode
    update_flumut=args.update_flumut_db
    off_apps=args.turn_off
    if args.single_sample.lower()=='on':
        single=True
    else:
        single=False
    if not single:
        flags['Master']['single']=False
    #### LOAD CONFIG PARAMETERS
    config=configparser.ConfigParser()
    config.read(config_file)
    if mode=='contig':
        print('Running AFluID in contig mode')
        flags['Master']['flumut']=False
        flags['Master']['genin']=False
        flags['Master']['nextclade']=False
        flags['Master']['getref']=True
    elif mode=='consensus':
        print('Running AFluID in consensus mode')
        flags['Master']['flumut']=True
        flags['Master']['genin']=True
        flags['Master']['nextclade']=True
        flags['Master']['getref']=False
    if off_apps:
        for i in off_apps:
            print(f'Turning off: {i}')
            flags['Master'][i]=False
#### SET PATHS
    samples=config['Paths']['samples']
    runs=config['Paths']['runs']
    references=config['Paths']['references']
    reports=config['Paths']['reports']
    logs=config['Paths']['logs']
    blasts=config['Paths']['blast_database']
    clusters=config['Paths']['cluster_database']
    metadata=config['Paths']['metadata']
    rm_previous=config['Functions']['remove_previous'] if\
          args.remove_previous else False
    cwd=os.getcwd()
    samples_p=os.path.join(cwd,samples)
    runs_p=os.path.join(cwd,runs)
    references_p=os.path.join(cwd,references)
    reports_p=os.path.join(cwd,reports)
    logs_p=os.path.join(cwd,logs)
    blasts_p=os.path.join(cwd,blasts)
    clusters_p=os.path.join(cwd,clusters)
    metadata_p=os.path.join(cwd,metadata)
    print('Single sample mode:',single)
    print('Remove previous files:',rm_previous)
    print('Configuration file:',config_file)
    print('Sample File:',filename)
    #### CHECK PATHS
    if os.path.exists(samples_p):
        if os.path.exists(runs_p):
            if os.path.exists(references_p):
                if os.path.exists(reports_p):
                    if os.path.exists(blasts_p):
                        if os.path.exists(clusters_p):
                            if os.path.exists(metadata_p):
                                print("Directories exist: Continuing analysis")
                            else:
                                print("Metadata path does not exist")
                                sys.exit()
                        else:
                            print("Cluster database path does not exist")
                            sys.exit()
                    else:
                        print("Blast database path does not exist")
                        sys.exit()
                else:
                    print("Reports path does not exist")
                    sys.exit()
            else:
                print("References path does not exist")
                sys.exit()
        else:
            print("Runs path does not exist")
            sys.exit()
    else:
        print("Samples path does not exist")
        sys.exit()
    if not os.path.exists(logs_p):
        os.makedirs(logs_p)
    if rm_previous:
        # Remove all files in RUN_DIR except *.pkl
        for file in glob.glob(os.path.join(runs_p, "*")):
            os.remove(file)


        # Remove specific files in SAMPLE_DIR
        for pattern in ["*_to_blast.fasta", "*_to_reblast.fasta", "*_to_flumut.fasta", "*_to_nextclade_H1.fasta",\
                         "*_to_nextclade_H3.fasta", "*_to_nextclade_H5.fasta","*_to_genin.fasta","format_*"]:
            for file in glob.glob(os.path.join(samples_p, pattern)):
                os.remove(file)

    ### LOAD PARAMETERS
    dict_cluster=pkl_load(os.path.join(clusters_p,config["Filenames"]["cluster_pkl"]))
    
    ###PIPELINE STEPS
    max=int(config["Sequence_Size"]["max"]) if not args.max_length else args.max_length
    min=int(config["Sequence_Size"]["min"]) if not args.min_length else args.min_length
    mappings=fasta_preprocess(filename,samples_p,min,max,flags,verbose=True)
    file=filename.replace('.fasta','')
    cd_hit_est_2d(os.path.join(samples_p,f'format_{filename}'),os.path.join(clusters_p,config["Filenames"]["cluster"]),os.path.join(runs_p,f'format_{file}'),float(config["CD-HIT"]["identity"]),logs_p,verbose=True)
    
    assignments=cluster_assign(f'format_{file}.clstr',os.path.join(samples_p,f'format_{filename}'),runs_p)
    cluster_compile(assignments,dict_cluster,os.path.join(runs_p,f'format_{file}.assign'),reports_p)
    cluster_report=cluster_miner(reports_p,file,samples_p,flags)
    if flags['Master']['C_BLAST']:
        blast_cluster=bclust(clusters_p,config["Filenames"]["cluster_metadata"],samples_p,blasts_p,config["Filenames"]["cluster"],runs_p,filename,flags,fasta=f"{filename.replace('.fasta','')}_to_blast.fasta")
        if flags['Master']['L_BLAST']:
            blast_report=reblast(metadata_p,config["Filenames"]["metadata"],samples_p,blasts_p,config["Filenames"]["l_blast"],runs_p,filename,fasta=f"{filename.replace('.fasta','')}_to_reblast.fasta")
            report_compiler(cluster_report,samples_p,filename,mappings,reports_p,flags,bclust_dict=blast_cluster,blast_dict=blast_report)
        else:
            report_compiler(cluster_report,samples_p,filename,mappings,reports_p,flags,bclust_dict=blast_cluster)
    else:
        report_compiler(cluster_report,samples_p,filename,mappings,reports_p,flags)
    #FORCING ANALYSIS
    
    if args.force:
        print('Forcing additional tools, please be mindful of the results')
        if 'flumut' in args.force:
           flags['Master']['flumut']=True
        if 'genin' in args.force:
            flags['Master']['genin']=True
        if 'getref' in args.force:  
            flags['Master']['getref']=True
    redirector(f"{file}_ID_Report.txt",flags,filename,mappings,samples_p,reports_p,force_flumut=flags['Master']['flumut'],force_genin=flags['Master']['genin'],force_getref=flags['Master']['getref'],mode=mode,single_sample=single)
    #flumut:
    if flags['Master']['flumut']:
        conform_to_flumut(flags,samples_p,filename)
        if update_flumut.upper()=='ON':
            update_flumut_db()
        run_flumut(reports_p,samples_p,file)
        remap_flumut_report(reports_p,mappings,file)
    #NextClade
    if flags['Master']['nextclade']:
        
        fasta_to_nextclade(flags,samples_p,filename)
        run_nextclade(f'{file}_to_nextclade',flags,reports_p,samples_p,file)
        if flags['Final Report']['Sequences for NextClade']['H1']!=[]:
            remap_nextclade(reports_p,f"{file}_H1_nextclade.tsv",mappings)
        if flags['Final Report']['Sequences for NextClade']['H3']!=[]:
            remap_nextclade(reports_p,f"{file}_H3_nextclade.tsv",mappings)
        if flags['Final Report']['Sequences for NextClade']['H5']!=[]:
            remap_nextclade(reports_p,f"{file}_H5_nextclade.tsv",mappings)
    
    #get reference:
    if flags['Master']['getref']:
        get_reference(flags['Final Report']['Get References'],references_p,config['Filenames']['ref_db'],filename,email=config['getref']['email'],update_local_db=True)
    
    #getting flagsdict into json
    with open(f'{os.path.join(reports_p,file)}_flags.json','w') as f:
        json.dump(flags,f)

if __name__=='__main__':
    main()
