#!/usr/bin/env python3
#### IMPORT DEPENDENCIES
import os
import sys
import argparse
import configparser
import subprocess
import json
import ast
import glob
import re
import numpy as np
from collections import defaultdict
from flu_utils import seq_get,dict_to_fasta,headers_from_mult_fas,parse_clstr,json_load,\
    seq_filter_get,mine_genotype_H,mine_genotype_N,convert_to_prop, mine_single_HA,mine_single_NA
from metadata_utils import ClusterReportTable, ClusterMetadata, SequenceMetadata
from structures import flagdict, taxa_dict, muts_loci_meaning , int_to_iupac, muts_interest,seg_lens, \
    iav_segments, contig_single_min_floor_bp, contig_single_min_step_bp
from final_report_utils import html_skeleton, generate_final_report, maybe_create_batch_artifacts
from copy import deepcopy
import pandas as pd
from gb_utils import fetch_genbank_list, append_genbank_from_list, get_gb, create_lookup,\
    write_gb,gb_to_fasta
from pathlib import Path
from datetime import datetime

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
    return f'>Seq-{counter}'

def fasta_preprocess(
    filename: str,
    samples_path: str,
    runs_path: str,
    file_tag: str,
    min_seq_len: int,
    max_seq_len: int,
    flags: dict,
    verbose: bool = False
) -> dict:
    """
    Preprocess a FASTA file by removing sequences outside of a length range and
    saving the processed sequences and header mappings.

    filename is the relative path inside samples_path.
    file_tag is the flat identifier used for generated files in runs_path.
    """
    fasta_dict = seq_get(os.path.join(samples_path, filename))

    to_pop = []
    for key in fasta_dict:
        if len(fasta_dict[key]) < min_seq_len or len(fasta_dict[key]) > max_seq_len:
            to_pop.append(key)

    for item in to_pop:
        if "Fasta Preprocess" in flags:
            flags["Fasta Preprocess"]["Rejected Sequences"].append(item)
        fasta_dict.pop(item)

    if verbose:
        if len(to_pop) == 1:
            print(f'{len(to_pop)} sequence was removed for being outside the length thresholds:')
        elif to_pop == []:
            print('All sequences were within length thresholds')
        else:
            print(f'{len(to_pop)} sequences were removed for being outside the length thresholds:')

    header_mapping = {}
    counter = 1
    for header in fasta_dict:
        key = seq_enum(counter)
        header_mapping[key] = header
        counter += 1

    mappings = {}
    for key in header_mapping:
        mappings[key] = fasta_dict[header_mapping[key]]

    outname_fasta = os.path.join(runs_path, f'format_{file_tag}')
    dict_to_fasta(mappings, outname_fasta)

    return header_mapping

#### CD-HIT
def cd_hit_est_2d(filename: str, cluster_reps: str, output: str, identity: float, logs_p: str, log_tag: str = None, verbose: bool = False) -> None:
    """
    Run CD-HIT-EST-2D to cluster sequences based on a given identity threshold.
    """
    tag = log_tag if log_tag is not None else os.path.basename(filename)

    if not verbose:
        log_fp = os.path.join(logs_p, f'{tag}_cd_hit.log')
        err_fp = os.path.join(logs_p, f'{tag}_cd_hit.err')
        with open(log_fp, 'a') as out, open(err_fp, 'a') as err:
            subprocess.run(
                ['cd-hit-est-2d', '-i', cluster_reps, '-i2', filename, '-o', output, '-c', str(identity), '-g', '1'],
                stdout=out,
                stderr=err
            )
    else:
        subprocess.run(['cd-hit-est-2d', '-i', cluster_reps, '-i2', filename, '-o', output, '-c', str(identity), '-g', '1'])

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
    with open(cl_assign, 'r') as fh:
        assign = fh.readlines()[1:]
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
    outpath=os.path.join(reports_path,f'{output_name.replace("format_","")}_clust_report.txt')
    
    #outputting report
    with open(f'{outpath}','w') as report:
        report.write('Sample_access\t%ID\tCluster\tCluster_rep\tGenotypes\tSegment\tHosts\tCountries\n')
        for sample in compiler:
            if len(compiler[sample])>2:
                report.write(f'{sample}\t{compiler[sample][0]}\t\
                            {compiler[sample][1]}\t{compiler[sample][6]}\t\
                            {compiler[sample][2]}\t{compiler[sample][3]}\t{compiler[sample][4]}\t{compiler[sample][5]}\n')
            elif len(compiler[sample])==0:
                continue  
            else:
                report.write(f'{sample}\t{compiler[sample][0]}\t{compiler[sample][1]}\n')
    print(f"Compiled report generated in {outpath}")

def cluster_miner(reports_p: str, file_tag: str, formatted_fasta_dir: str, runs_p: str, flags: dict) -> dict:
    """
    Mine cluster report for unassigned samples and output them to BLAST.

    Parameters
    ----------
    reports_p : str
        Path to the reports directory.
    file_tag : str
        Sample file base name.
    formatted_fasta_dir : str
        Path to the samples directory.
    runs_p : str
        Path to the runs directory.
    flags : dict
        Flags dictionary to update.

    Returns
    -------
    dict
        Dictionary mapping sequence identifiers to the same 8 output fields
        already used in the original function.
    """
    report_fp = os.path.join(reports_p, f"{file_tag}_clust_report.txt")
    fasta_fp = os.path.join(formatted_fasta_dir, f"format_{file_tag}.fasta")
    blast_fp = os.path.join(runs_p, f"{file_tag}_to_blast")
    # Load report object
    clust_rep = ClusterReportTable(report_fp)
    clust_rep.convert_to_prop("genotypes")
    clust_rep.convert_to_prop("hosts")
    clust_rep.convert_to_prop("countries")
    hs = clust_rep.extract_h()
    ns = clust_rep.extract_n()

    # Load tabular report
    evaluate = pd.read_table(report_fp, index_col=False)

    if "Sample_access" not in evaluate.columns or "Cluster" not in evaluate.columns:
        raise ValueError("The cluster report must contain 'Sample_access' and 'Cluster' columns.")

    # Safe conversion to plain Python list
    to_blast = (
        evaluate.loc[evaluate["Cluster"] == ">Unassigned", "Sample_access"]
        .dropna()
        .astype(str)
        .tolist()
    )
    to_blast_set = set(to_blast)

    to_report = {}

    for key, row in clust_rep.data.items():
        key_str = str(key)

        # Skip unassigned sequences
        if key_str in to_blast_set:
            continue

        # Protect against short/incomplete rows without changing output shape
        safe_row = list(row) if row is not None else []
        safe_row += [""] * max(0, 7 - len(safe_row))

        to_report[key_str] = [
            str(safe_row[2]).replace(" ", "") if safe_row[2] is not None else "",
            str(safe_row[1]).replace("                        >", ">") if safe_row[1] is not None else "",
            str(safe_row[0]).replace("%", "") if safe_row[0] is not None else "",
            safe_row[4] if safe_row[4] is not None else "",
            safe_row[3] if safe_row[3] is not None else "",
            safe_row[5] if safe_row[5] is not None else "",
            safe_row[6] if safe_row[6] is not None else "",
            "CD-HIT"
        ]

    # Update flags and export FASTA for BLAST
    flags.setdefault("CD-HIT", {})
    flags["CD-HIT"].setdefault("Unclustered Sequences", [])
    flags.setdefault("Master", {})

    if to_blast:
        flags["CD-HIT"]["Unclustered Sequences"].extend(to_blast)
        to_blast_dict = seq_filter_get(fasta_fp, to_blast)
        dict_to_fasta(to_blast_dict, blast_fp)
        flags["Master"]["C_BLAST"] = True
    else:
        flags["Master"]["C_BLAST"] = False

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

def bclust(metadata_p:str,metadata_f:str,samples_p:str,blast_p:str,blast_db:str,runs_p:str,file_tag:str,flags:dict,num_threads:int=2,max_tar_seq:int=7,fasta:str='to_blast.fasta') -> dict:
    '''
    Run BLAST analysis on unassigned samples.
    Parameters:
    metadata_p (str): The path to the metadata directory.
    samples_p (str): The path to the samples directory.
    blast_p (str): The path to the BLAST database directory.
    runs_p (str): The path to the runs directory.
    filename (str): The name of the sample file.
    num_threads (int): The number of threads to use for BLAST analysis.
    max_tar_seq (int): The maximum number of target sequences to return.
    fasta (str): The name of the FASTA file to be analyzed
    Outputs a .fasta file for BLAST analysis.
    '''
    access=headers_from_mult_fas([os.path.join(samples_p,fasta)],only_name=True)
    queries=[key for key in access]
    subprocess.run(['blastn', '-db', os.path.join(blast_p,blast_db), '-query', os.path.join(samples_p,fasta),
               '-out', os.path.join(runs_p,f"{file_tag}_bcrun.txt"), '-outfmt', '6', '-num_threads', str(num_threads), '-max_target_seqs', str(max_tar_seq)])
    b_tab=best_blast(runs_p,f"{file_tag}_bcrun.txt",thres_id=90.0)
    pd.set_option('display.max_colwidth', None)
    metadata=pd.read_table(os.path.join(metadata_p,metadata_f),index_col=False)
    assigned=set(b_tab.keys())
    queries=set(queries)
    to_reblast=queries-assigned
    for key in b_tab:
        b_tab[key].append(metadata['Cluster'][metadata['Representative']==b_tab[key][0]].to_string(index=False,min_rows=None,max_rows=None))
        b_tab[key].append(metadata['Segments'][metadata['Representative']==b_tab[key][0]].to_string(index=False,min_rows=None,max_rows=None))
        b_tab[key][3]=ast.literal_eval(b_tab[key][3])
        b_tab[key].append(metadata['Genotypes'][metadata['Representative']==b_tab[key][0]].to_string(index=False,min_rows=None,max_rows=None))
        b_tab[key][4]=ast.literal_eval(b_tab[key][4])
        b_tab[key].append(metadata['Hosts'][metadata['Representative']==b_tab[key][0]].to_string(index=False,min_rows=None,max_rows=None))
        b_tab[key][5]=ast.literal_eval(b_tab[key][5])
        b_tab[key].append(metadata['Countries'][metadata['Representative']==b_tab[key][0]].to_string(index=False,min_rows=None,max_rows=None))
        b_tab[key][6]=ast.literal_eval(b_tab[key][6])

    if list(to_reblast) != []:
        for i in to_reblast:
            flags['BLAST']['Sequences unassigned against cluster representatives'].append(i)
        to_reblast_dict=seq_filter_get(os.path.join(samples_p,fasta),to_reblast)
        dict_to_fasta(to_reblast_dict,os.path.join(runs_p,f'{file_tag}_to_reblast'))
        flags['Master']['L_BLAST']=True
    else:
        flags['Master']['L_BLAST']=False
    return b_tab

def reblast(metadata_p:str,metadata_f:str,samples_p:str,blast_p:str,blast_db:str,runs_p:str,file_tag:str,threads:int=2,max_tar_seq:int=7,fasta:str='to_reblast.fasta') -> dict:

    '''
    Re-run BLAST analysis on unassigned samples.
    Parameters:
    metadata_p (str): The path to the metadata directory.
    samples_p (str): The path to the samples directory.
    blast_p (str): The path to the BLAST database directory.
    runs_p (str): The path to the runs directory.
    filename (str): The name of the sample file.
    num_threads (int): The number of threads to use for BLAST analysis.
    max_tar_seq (int): The maximum number of target sequences to return.
    fasta (str): The name of the FASTA file to be analyzed

    '''
    subprocess.run(['blastn', '-db', os.path.join(blast_p,blast_db), '-query', os.path.join(samples_p,fasta),
               '-out', os.path.join(runs_p,f"{file_tag}_brun.txt"), '-outfmt', '6', '-num_threads', str(threads), '-max_target_seqs', str(max_tar_seq)])
    report=best_blast(runs_p,f"{file_tag}_brun.txt")
    metadata=pd.read_csv(os.path.join(metadata_p,metadata_f),sep=';',index_col=False)
    for key in report:
        report[key].append(metadata['SEGMENT'][metadata['ACCESSION']==report[key][0]].to_string(index=False))
        report[key].append(metadata['GENOTYPE'][metadata['ACCESSION']==report[key][0]].to_string(index=False))
        report[key].append(metadata['HOST'][metadata['ACCESSION']==report[key][0]].to_string(index=False))
        report[key].append(metadata['COUNTRY'][metadata['ACCESSION']==report[key][0]].to_string(index=False))
    return report

#### REPORT COMPILER
def report_compiler(clust_dict, samples_p, file_tag, mappings_dict, reports_p, flags, bclust_dict={}, blast_dict={}):
    '''
    Compile BLAST and cluster reports into a final report.
    Parameters:
    clust_dict (dict): The cluster report dictionary.
    samples_p (str): Path to samples file
    bclust_dict (dict): The BLAST cluster report dictionary.
    blast_dict (dict): The BLAST report dictionary.
    file_tag (str): The name of the sample file.
    mappings_dict (dict): The mappings dictionary.
    reports_p (str): The path to the reports directory.
    Outputs a final report file.
    '''
    seqs = headers_from_mult_fas(
    [os.path.join(samples_p, f'format_{file_tag}.fasta')],
    only_name=True,
    out_list=True)
    to_report=clust_dict
    if bclust_dict:
        for key in bclust_dict:
            if key not in to_report:
                to_report[key]=[]
                to_report[key].append(bclust_dict[key][0])
                to_report[key].append(bclust_dict[key][2])
                to_report[key].append(bclust_dict[key][1])
                to_report[key].append(bclust_dict[key][3])
                to_report[key].append(convert_to_prop(bclust_dict[key][4]))
                to_report[key].append(convert_to_prop(bclust_dict[key][5]))
                to_report[key].append(convert_to_prop(bclust_dict[key][6]))
                to_report[key].append('C-BLAST')
    if blast_dict:
        for key in blast_dict:
            if key not in to_report:
                to_report[key]=[]
                to_report[key].append(blast_dict[key][0])
                to_report[key].append('UNCLUSTERED')
                to_report[key].append(blast_dict[key][1])
                to_report[key].append(blast_dict[key][2].replace(' ',''))
                to_report[key].append(blast_dict[key][3].replace(' ',''))
                to_report[key].append(blast_dict[key][4])
                to_report[key].append(blast_dict[key][5].replace(' ',''))
                to_report[key].append('L-BLAST')
    assigned=set(to_report.keys())
    seqs=set(seqs)
    to_remote=seqs-assigned
    mappings=mappings_dict
    if list(to_remote)!=[]:
        for seq in to_remote:
            to_report[seq]=['Unassigned','NA','NA','NA','NA','NA','NA','NA']
    with open(os.path.join(reports_p, f'{file_tag}_ID_Report.txt'), 'w') as report:
        report.write('Sample_name\tRepresentative\tCluster\t%ID\tSegment\tGenotype\tHost\tCountry\tAssigned_by\n'.upper())
        for key in to_report:
            mapped=mappings[f'>{key}']
            #print(to_report[key])
            report.write(f'{mapped}\t{to_report[key][0]}\t\
                         {to_report[key][1]}\t{to_report[key][2]}\t{to_report[key][3]}\t\
                            {to_report[key][4]}\t{to_report[key][5]}\t{to_report[key][6]}\t{to_report[key][7]}\n')
    if len(list(to_remote))>0:
        flags['BLAST']['Sequences unassigned against local database'].extend(to_remote)
    #print(f'Report generated in {os.path.join(reports_p,filename.replace(".fasta",""))}_ID_Report.txt')

#### REPORT PARSERS
def _parse_segment_cell(x):
    """Parse SEGMENT cell from ID_Report into int 1..8 (or None)."""
    if pd.isna(x):
        return None
    s = str(x).strip()
    if not s:
        return None
    try:
        v = ast.literal_eval(s) if s.startswith('{') else s
    except Exception:
        v = s
    if isinstance(v, dict) and v:
        k = next(iter(v.keys()))
        try:
            return int(k)
        except Exception:
            return None
    try:
        return int(v)
    except Exception:
        return None


def _parse_pident_cell(x):
    """Parse %ID cell from ID_Report into float (or None)."""
    if pd.isna(x):
        return None
    s = str(x).strip().replace('%', '')
    if not s or s.upper() == 'NA':
        return None
    try:
        return float(s)
    except Exception:
        return None


def compute_best_segments_iterative(id_report_fp: str, formatted_fasta_fp: str, mappings: dict, start_min_len: int,
                                   step_bp: int = contig_single_min_step_bp,
                                   floor_bp: int = contig_single_min_floor_bp) -> tuple:
    """Select a single best contig per segment without overwriting.

    Identification is assumed to have been run once on the full formatted FASTA.
    We then iterate selection thresholds: start_min_len, start_min_len-step, ... until floor_bp.

    Selection order: contig length desc, then %ID desc.
    Fill-once: once a segment is filled, it is never replaced.

    Returns:
      (best_segments_internal, min_len_used)
    """
    best = {seg: None for seg in iav_segments}

    if not os.path.exists(id_report_fp) or not os.path.exists(formatted_fasta_fp):
        return best, int(floor_bp)

    df = pd.read_table(id_report_fp, index_col=False)
    if df is None or df.empty:
        return best, int(floor_bp)

    if 'SAMPLE_NAME' not in df.columns or 'SEGMENT' not in df.columns:
        return best, int(floor_bp)

    # Original header -> internal id (e.g. ">Seq_1")
    remap = {orig: internal for internal, orig in mappings.items()}
    seqs = seq_get(formatted_fasta_fp)

    pident_col = '%ID' if '%ID' in df.columns else None
    df['segment_num'] = df['SEGMENT'].apply(_parse_segment_cell)
    df['pident'] = df[pident_col].apply(_parse_pident_cell) if pident_col else None
    df['pident'] = df['pident'].fillna(-1.0).astype(float)
    df['internal_id'] = df['SAMPLE_NAME'].astype(str).map(remap)

    def _len_for_internal(sid):
        if not isinstance(sid, str):
            return 0
        seq = seqs.get(sid)
        return len(seq) if seq is not None else 0

    df['contig_len'] = df['internal_id'].apply(_len_for_internal)

    cand = df[(df['segment_num'].isin(range(1, 9))) & (df['contig_len'] > 0) & (df['internal_id'].notna())].copy()
    if cand.empty:
        return best, int(floor_bp)

    cand = cand.sort_values(['contig_len', 'pident'], ascending=[False, False])
    candidates = list(cand.itertuples(index=False))

    start = int(start_min_len) if start_min_len is not None else int(floor_bp)
    start = max(start, int(floor_bp))
    step = max(1, int(step_bp))
    floor = int(floor_bp)

    thresholds = [start] if start == floor else list(range(start, floor, -step)) + [floor]
    min_used = thresholds[-1]

    for th in thresholds:
        min_used = th
        for row in candidates:
            if int(row.contig_len) < th:
                continue
            seg_num = int(row.segment_num) if row.segment_num is not None else None
            seg_name = int_to_iupac.get(seg_num)
            if not seg_name or seg_name not in best:
                continue
            if best[seg_name] is None:
                best[seg_name] = str(row.internal_id)
        if all(best.values()):
            break

    return best, int(min_used)

#### REDIRECTOR
def safe_literal_eval(x):
    if pd.isna(x):
        return x
    s = str(x).strip()
    if not s:
        return x
    try:
        return ast.literal_eval(s)
    except Exception:
        return x
def redirector(report, flags, file_tag, mappings, runs_p, reports_p,
               force_flumut, force_genin, force_getref, mode, single_sample=True) -> None:
    '''
    Redirects samples to aditional post-identification tools.
    Linchpin function that populates the 'Sample' and 'Final_report' dicts of flagsdict
    Parameters:
    report (str): The path to the final report file.
    flags (dict): The flags dictionary.

    '''
    #opening report and extracting data
    report_df=pd.read_table(os.path.join(reports_p,report), index_col=False)
# contig + single-sample: restrict downstream work to the selected best contigs
    if mode == 'contig' and single_sample:
        best = flags.get('Sample', {}).get('best_segments', {})

        if isinstance(best, dict):
            keep_internal = {v for v in best.values() if v}
            keep_original = {mappings[iid] for iid in keep_internal if iid in mappings}

            #print("keep_internal:", keep_internal)
            print("keep_original:", keep_original)

        if 'SAMPLE_NAME' in report_df.columns:
            if keep_original:
                report_df = report_df.loc[report_df['SAMPLE_NAME'].isin(keep_original)].copy()
                print(f"Filtered report to {len(report_df)} rows")
            else:
                print("keep_original empty: keeping original report unchanged")
        report_df.to_csv(os.path.join(reports_p,report),sep='\t',index=False)
    report_df=report_df.dropna()
    report_df['SEGMENT'] = report_df['SEGMENT'].apply(safe_literal_eval)
    report_df['SEGMENT']=report_df['SEGMENT'].apply(lambda x: list(x.keys())[0] if type(x)==dict else int(x))
    Segments=report_df['SEGMENT'].to_list()
    Genotypes=report_df['GENOTYPE'].to_list()
    Seq_names=report_df['SAMPLE_NAME'].to_list()
    #Type conversions and dictionary creations
    for i in range(len(Genotypes)):
        try:
            Genotypes[i] = safe_literal_eval(Genotypes[i])
        except NameError:
            continue
        Genotypes[i]=list(Genotypes[i].keys()) if\
            type(Genotypes[i])==dict else Genotypes[i]
    Segments=[int(i) for i in Segments]
    Seq_gen=dict(zip(Seq_names,Genotypes))
    Seq_seg=dict(zip(Seq_names,Segments))
    Seq_ref = dict(zip(report_df['SAMPLE_NAME'], report_df['REPRESENTATIVE']))
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
            force_flumut = True
        if flags['Sample']['H_gen']=='H5' or force_flumut:
            for i in Seq_names:
                flags['Final Report']['Sequences for FluMut'].append(remap[i])

    #Opening formatted fasta 
        if single_sample:
            formatted = seq_get(os.path.join(runs_p, f'format_{file_tag}.fasta'))
            for seg in ('HA','NA','PB2','PB1','PA','NP','NS','MP'):
                try:
                    flags['Sample'][f'{seg}_len']=len(formatted[flags['Sample'][seg][0]])
                except IndexError:
                    flags['Sample'][f'{seg}_len']=0
                except KeyError:
                    flags['Sample'][f'{seg}_len']=0
                    
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
    subprocess.run(['flumut', '--update'])

def conform_to_flumut(flagdict: dict, sample_path: str, file_tag: str) -> None:
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
    flt_fasta = seq_get(os.path.join(sample_path, f'format_{file_tag}.fasta'))
    outdict = {}
    for key in flt_fasta:
        if key in seq_seg:
            outdict[f'{key}_{seq_seg[key]}'] = flt_fasta[key]
    dict_to_fasta(outdict, os.path.join(sample_path, f'{file_tag}_to_flumut'))
        
def run_flumut(reports_p: str, samples_p: str, file_tag: str, regex="(.+)\_(.+)"):
    """
    Executes the flumut tool to analyze sequence data and generate mutation reports.

    Parameters:
    reports_p (str): The path to the directory where the report files will be stored.
    samples_p (str): The path to the directory containing the sample FASTA file.
    file_tag (str): The base name of the sample file, used to generate report filenames.
    regex (str): A regular expression pattern used by flumut for parsing sequence headers. Defaults to "(.+)\|(.+)".

    Returns:
    None: This function does not return a value. It executes a system command to run flumut and generate output files.
    """
    print("Running flumut")
    subprocess.run(['flumut', '-m', os.path.join(reports_p,f"{file_tag}_markers.tsv"), '-M', os.path.join(reports_p,f"{file_tag}_mutations.tsv"),
                    '-l', os.path.join(reports_p,f"{file_tag}_literature.tsv"), os.path.join(samples_p,f"{file_tag}_to_flumut.fasta"), '-n', regex])
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
    #Create a mapping dictionary with the original sequence names replaced by the new names
    mapping_dict = {key.replace('>', ''): value for key, value in mappings.items()}

    #Read the markers report file and replace the original sequence names with the new names
    mk_report = pd.read_table(os.path.join(reports_dir,f'{filename}_markers.tsv'))
    mk_report = mk_report.replace(mapping_dict)

    #Read the mutations report file and replace the original sequence names with the new names
    mut_report = pd.read_table(os.path.join(reports_dir,f'{filename}_mutations.tsv'), low_memory=False)
    mut_report = mut_report.replace(mapping_dict)

    #Save the modified markers report file9
    mk_report.to_csv(os.path.join(reports_dir,f'{filename}_markers.tsv'), sep='\t', index=False)

    #Save the modified mutations report file
    mut_report.to_csv(os.path.join(reports_dir,f'{filename}_mutations.tsv'), sep='\t', index=False)
def mut_miner(dataframe,muts_of_interest,flagsdict):
    """
    Mines the flumut markers report for mutations of interest according to Alvarez et al 2025 and 
    Mohapra et al 2023.
    Parameters:
    dataframe (str): The path to the flumut markers report file.
    muts_of_interest (dict): A dictionary containing mutations of interest for each segment.
    flagsdict (dict): A dictionary to store the mined mutations of interest.
    Returns:
    None: This function does not return a value. It updates the flagsdict with the mined 
    """
    # Load mutations of interest from pickle file
    # MUTATION PATTERN
    mut_patt=r'([A-Z])(\d+)([A-Z])'
    # Read the DataFrame
    df=pd.read_table(dataframe, index_col=False)
    samples=list(df['Sample'].unique())
    sample_muts={}
    for i in samples:
        sample_muts[i]=set()
        list_i=df['Mutations in your sample'][df['Sample']==i].to_list()
        for j in list_i:
            sample_muts[i].add(j)
        sample_muts[i]=list(sample_muts[i])
    sample_seg={sample:set() for sample in sample_muts}
    for sample in sample_muts:
        for i in range(len(sample_muts[sample])):
            sample_seg[sample].add(sample_muts[sample][i].split(':')[0])
            sample_muts[sample][i]=sample_muts[sample][i].split(':')[1]
    loci_seg={'PB1-F2': 'PB1', 'PA-X': 'PA', 'HA1-5': 'HA', 'HA2-5': 'HA', 'NA-1':'NA', 'NA-2':'NA','M1':'MP', 'M2':'MP', 'NS-1':'NS', 'NS-2':'NS'}
    for sample in sample_seg:
        sample_seg[sample]=list(sample_seg[sample])
    for sample in sample_seg:
        for i in sample_seg[sample]:
            if i in loci_seg.keys():
                sample_seg[sample].append(loci_seg[i])
                sample_seg[sample].remove(i)
    for i in sample_seg:
        for j in sample_seg[i]:
            if j in loci_seg.keys():
                sample_seg[i].remove(j)
    for sample in sample_seg:
        sample_seg[sample]=set(sample_seg[sample])
    for sample in sample_muts:
        for i in range(len(sample_muts[sample])):
            ex=re.match(mut_patt,sample_muts[sample][i])
            if ex:
                sample_muts[sample][i]=f'{ex.group(2)}{ex.group(3)}'
    sample_muts_of_interest={sample: {seg: set() for seg in sample_seg[sample]} for sample in sample_muts}
    for sample in sample_muts_of_interest:
        for segment in sample_muts_of_interest[sample]:
            if segment in muts_of_interest.keys():    
                for mut in sample_muts[sample]:
                    if mut in muts_of_interest[segment]:
                        sample_muts_of_interest[sample][segment].add(mut)
    for i in sample_muts_of_interest:
        for j in sample_muts_of_interest[i]:
            sample_muts_of_interest[i][j]=list(sample_muts_of_interest[i][j])
            #flagsdict["Sample"][f"{j}_muts"]=[]
            for mut in sample_muts_of_interest[i][j]:
                #print(f'Sample: {i} | Segment: {j} | Mutation: {mut}')
                flagsdict["Sample"][f"{j}_muts"].append(mut)

#### GET REFERENCE PATH
def get_reference(listref,references_p,db_p,file_tag,email=None, update_local_db=False):
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
    outname = f'{file_tag}_references'
    write_gb(output,os.path.join(references_p,f'{outname}.gb'))
    gb_to_fasta(os.path.join(references_p,f'{outname}.gb'),os.path.join(references_p,f'{outname}.fasta'))

    if update_local_db:
        append_genbank_from_list(os.path.join(references_p,db_p),remote)
#### NEXTCLADE PATH
def fasta_to_nextclade(flagsdict, samples_p, file_tag):
    to_nextclade = {'H1': [], 'H3': [], 'H5': []}

    for genotype in to_nextclade:
        to_nextclade[genotype].extend(
            flagsdict['Final Report']['Sequences for NextClade'][genotype]
        )

    flt_fasta = seq_get(os.path.join(samples_p, f'format_{file_tag}.fasta'))

    for genotype, seqs in to_nextclade.items():
        if not seqs:
            continue

        outdict = {}
        for sample in seqs:
            if sample in flt_fasta:
                outdict[sample] = flt_fasta[sample]

        if outdict:
            dict_to_fasta(outdict, os.path.join(samples_p, f'{file_tag}_to_nextclade_{genotype}'))
def run_nextclade(fasta,flagsdict,reports_p,samples_p,file_tag):
    dict_builds_broad={'H1':'flu_h1n1pdm_ha_broad','H3':'flu_h3n2_ha_broad','H5':'community/moncla-lab/iav-h5/ha/all-clades'}
    for key in flagsdict['Final Report']['Sequences for NextClade']:
        if flagsdict['Final Report']['Sequences for NextClade'][key] != []:
            fasta_gen=key    
            subprocess.run(['nextclade', 'run', '-d', dict_builds_broad[fasta_gen], '--output-tsv', os.path.join(reports_p,f"{file_tag}_{fasta_gen}_nextclade.tsv"), os.path.join(samples_p,f"{fasta}_{fasta_gen}.fasta")])
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

#### GENIN PATH

def mine_clade(clade_file, flagsdict, mappings_dict):
    """
    Mines the clade from the NextClade report file based on the HA genotype.
    """
    clade = ""

    if (
        flagsdict["Sample"]["H_gen"] in ["H1", "H3", "H5"]
        and len(set(flagsdict["Final Report"]["Sequences for NextClade"][flagsdict["Sample"]["H_gen"]])) == 1
    ):
        seq = mappings_dict[flagsdict['Sample']['HA'][-1]]
        nc = pd.read_table(clade_file, index_col=False)
        nc["seqName"] = nc["seqName"].astype(str).str.strip()

        sel = nc.loc[nc["seqName"] == seq]

        if sel.empty:
            clade = "Unassigned"
        else:
            error_val = sel["errors"].iloc[0] if "errors" in sel.columns else None
            clade_val = sel["clade"].iloc[0] if "clade" in sel.columns else None

            if pd.notna(error_val) and str(error_val).strip() != "":
                clade = "Unassigned by error"
            elif pd.isna(clade_val) or str(clade_val).strip() == "":
                clade = "Unassigned"
            else:
                clade = str(clade_val).strip()
    else:
        clade = "Non-applicable"

    return clade

def to_genin2(flagsdict):
    """
    Determines if the sample should be directed to GenIn2 based on HA genotype and clade
    Parameters:
    flagsdict (dict): A dictionary containing flags and sequence information, including HA genotype and clade.
    Returns:
    bool: True if the sample should be directed to GenIn2, False otherwise."""
    if flagsdict["Sample"]["H_gen"] == "H5":
        if flagsdict["Sample"]["clade"] == "2.3.4.4b":
            return True
    return False

def run_genin2(samples_p, file_tag, genin_report):
    """
    Executes the GenIn2 tool for genotyping analysis based on HA genotype and clade.
    Parameters:
    flagsdict (dict): A dictionary containing flags and sequence information, including HA genotype and clade.
    samples_p (str): The path to the directory containing the sample FASTA file.
    filename (str): The name of the FASTA file to be processed
    genin_report (str): The path to the directory where the report files will be stored.
    Returns:
    None: This function does not return a value. It executes a system command to run GenIn2 and generate output files.
    """
    subprocess.run(['genin2','-o', genin_report, os.path.join(samples_p,f"{file_tag}_to_flumut.fasta")])

def remap_genin2(reports_p, genin_report, mappings, flagsdict):
    genin_path = os.path.join(reports_p, genin_report)
    genin_df = pd.read_table(genin_path, index_col=False)

    genin_df['Sample Name'] = genin_df['Sample Name'].apply(
        lambda x: str(mappings[f'>{x}']).strip() if f'>{x}' in mappings else str(x).strip()
    )

    genin_df.to_csv(genin_path, sep='\t', index=False)

    segment_cols = ['PB2', 'PB1', 'PA', 'NP', 'NA', 'MP', 'NS']
    result = {}
    skipped = []

    for _, row in genin_df.iterrows():
        sample_name = str(row.get('Sample Name', '')).strip()
        genotype = str(row.get('Genotype', '')).strip()

        sub_genotype = row.get('Sub-genotype', '')
        if pd.isna(sub_genotype) or str(sub_genotype).strip() == '':
            sub_genotype = 'N/A'
        else:
            sub_genotype = str(sub_genotype).strip()

        valid = row[segment_cols][row[segment_cols].notna() & (row[segment_cols] != '?')]

        if len(valid) == 0:
            skipped.append({'sample': sample_name, 'reason': 'no_valid_segment'})
            continue

        if len(valid) > 1:
            skipped.append({'sample': sample_name, 'reason': 'multiple_valid_segments', 'values': valid.to_dict()})
            continue

        segment_name = str(valid.index[0]).strip()
        segment_value = str(valid.iloc[0]).strip()

        result[sample_name] = {
            'segment': segment_name,
            'Genotype': genotype,
            'Sub-genotype': sub_genotype,
            'value': segment_value
        }

    flagsdict['Sample']['Genin_genotypes'] = result
    flagsdict['Sample']['Genin_skipped_rows'] = skipped
    return flagsdict


#### MAIN
def parser():
    parser = argparse.ArgumentParser(
        prog='AFluID',
        description='AFluID: Automated Influenza Identification Pipeline'
    )
    parser.add_argument('-c', '--config', type=str, help='Configuration file path', default=os.path.join(os.path.dirname(__file__), 'config.ini'), required=False)
    parser.add_argument('-f', '--filename', type=str, help='FASTA file name', required=False)
    parser.add_argument('-m', '--mode', type=str, help='Mode of operation', choices=('contig', 'consensus'), required=True)
    parser.add_argument('-b', '--batch', action='store_true', help='Run pipeline for all FASTA files in batch directory')
    parser.add_argument('-bd', '--batch_dir', type=str, default='', help='Subdirectory inside samples/ to use for batch mode')
    parser.add_argument('-ff', '--force', help='Force run of aditional tools', nargs='*', choices=('flumut', 'genin', 'getref'))
    parser.add_argument('-fdb', '--update_flumut_db', type=str, help='turn off auto-update for flumut db', choices=('on', 'off'), default='on', required=False)
    parser.add_argument('-ss', '--single_sample', type=str, help='Single sample mode', default='on', choices=('on', 'off'), required=False)
    parser.add_argument('-off', '--turn_off', help='Turn Off additional analysis tools', nargs='*', choices=('flumut', 'genin', 'nextclade', 'getref'))
    parser.add_argument('-rm', '--remove_previous', type=str, help='Remove previous files', default='on', choices=('on', 'off'), required=False)
    parser.add_argument('-Ml', '--max_length', type=int, help='Maximum sequence length', required=False)
    parser.add_argument('-ml', '--min_length', type=int, help='Minimum sequence length', required=False)

    args = parser.parse_args()

    if not args.batch and not args.filename:
        parser.error("--filename is required unless --batch is used")

    return args
def make_sample_stem(filename: str) -> str:
    """
    Get only the FASTA basename without extension.

    Examples
    --------
    H5N1_IRMA/PT_A_H5N1_SAMPLE_01_IRMA.fasta -> PT_A_H5N1_SAMPLE_01_IRMA
    sample1.fasta -> sample1
    """
    return Path(filename).stem


def make_run_prefix(batch_dir: str = '', timestamp: str = '') -> str:
    """
    Build the batch run prefix.

    Examples
    --------
    batch_dir='H5N1_IRMA', timestamp='20260420_161518'
    -> H5N1_IRMA_20260420_161518
    """
    batch_label = batch_dir.strip().replace(os.sep, "__").replace("/", "__") if batch_dir else "batch_root"
    return f"{batch_label}_{timestamp}" if timestamp else batch_label


def make_output_stem(filename: str, batch: bool = False, batch_dir: str = '', timestamp: str = '') -> str:
    """
    Build the stem used for generated output files.

    Single-file mode:
        PT_A_H5N1_SAMPLE_01_IRMA

    Batch mode:
        H5N1_IRMA_20260420_161518_PT_A_H5N1_SAMPLE_01_IRMA
    """
    sample_stem = make_sample_stem(filename)

    if batch:
        run_prefix = make_run_prefix(batch_dir=batch_dir, timestamp=timestamp)
        return f"{run_prefix}_{sample_stem}"

    return sample_stem


def make_file_tag(filename: str) -> str:
    """
    Convert a relative FASTA path into a safe flat tag for output filenames.

    Examples
    --------
    sample1.fasta -> sample1
    projectA/run3/sample2.fasta -> projectA__run3__sample2
    """
    p = Path(filename)
    no_suffix = p.with_suffix("")
    return str(no_suffix).replace(os.sep, "__")


def find_fasta_files(samples_p: str, batch_dir: str = '', recursive: bool = True):
    """
    Discover FASTA files inside samples_p or inside a given subdirectory of samples_p.

    Parameters
    ----------
    samples_p : str
        Root samples directory from config.
    batch_dir : str
        Optional subdirectory inside samples_p.
    recursive : bool
        Whether to search recursively.

    Returns
    -------
    list[str]
        Relative paths from samples_p.
    """
    samples_root = Path(samples_p)

    if batch_dir:
        search_root = samples_root / batch_dir
    else:
        search_root = samples_root

    if not search_root.exists():
        raise FileNotFoundError(f"Batch directory does not exist: {search_root}")

    if not search_root.is_dir():
        raise NotADirectoryError(f"Batch path is not a directory: {search_root}")

    patterns = ["*.fasta", "*.fa", "*.fna"]
    found = []

    if recursive:
        for pattern in patterns:
            found.extend(search_root.rglob(pattern))
    else:
        for pattern in patterns:
            found.extend(search_root.glob(pattern))

    relpaths = [str(p.relative_to(samples_root)) for p in found if p.is_file()]
    return sorted(set(relpaths))

def run_pipeline_for_file(
    filename: str,
    flags_template: dict,
    config,
    samples_p: str,
    runs_p: str,
    references_p: str,
    reports_p: str,
    logs_p: str,
    blasts_p: str,
    clusters_p: str,
    metadata_p: str,
    mode: str,
    single: bool,
    update_flumut: str,
    off_apps,
    force_apps,
    rm_previous: bool,
    max_len: int,
    min_len: int,
    batch: bool = False,
    batch_dir: str = '',
    timestamp: str = ''
):
    """
    Run full AFluID pipeline for one FASTA file.

    Parameters
    ----------
    filename : str
        Relative path from samples_p to the FASTA file.
    batch : bool, default False
        Whether this file is being processed under batch mode.
    batch_dir : str, default ''
        Batch subdirectory label to include in output names.
    timestamp : str, default ''
        Run timestamp to include in batch output names.
    """
    flags = deepcopy(flags_template)
    muts_of_interest = muts_interest
    flags['Master']['single'] = single


    if mode == 'contig':
        print('Running AFluID in contig mode')
        flags['Master']['flumut'] = False
        flags['Master']['genin'] = False
        flags['Master']['nextclade'] = False
        flags['Master']['getref'] = True
    elif mode == 'consensus':
        print('Running AFluID in consensus mode')
        flags['Master']['flumut'] = True
        flags['Master']['genin'] = False
        flags['Master']['nextclade'] = True
        flags['Master']['getref'] = False

    if off_apps:
        for i in off_apps:
            print(f'Turning off: {i}')
            flags['Master'][i] = False
    forced_genin=False
    if force_apps:
        print('Forcing additional tools, please be mindful of the results')
        forced_flumut = force_apps is not None and 'flumut' in force_apps
        forced_genin = force_apps is not None and 'genin' in force_apps
        forced_getref = force_apps is not None and 'getref' in force_apps
        if force_apps:
            print('Forcing additional tools, please be mindful of the results')
            if forced_flumut:
                flags['Master']['flumut'] = True
            if forced_genin:
                flags['Master']['genin'] = True
            if forced_getref:
                flags['Master']['getref'] = True

    selection_start_min = min_len

    print('-' * 80)
    print('Processing file:', filename)
    print('Single sample mode:', single)
    print('Max sequence length:', max_len)
    print('Min sequence length:', min_len)

    sample_stem = Path(filename).stem

    if batch:
        batch_label = batch_dir.strip().replace(os.sep, "_").replace("/", "_") if batch_dir else "batch"
        if timestamp:
            output_tag = f"{batch_label}_{timestamp}_{sample_stem}"
        else:
            output_tag = f"{batch_label}_{sample_stem}"
    else:
        output_tag = sample_stem

    print('Output tag:', output_tag)

    formatted_fasta_base = os.path.join(runs_p, f'format_{output_tag}')
    formatted_fasta_fp = f'{formatted_fasta_base}.fasta'

    preprocess_min = min_len
    if mode == 'contig' and single:
        preprocess_min = contig_single_min_floor_bp

    mappings = fasta_preprocess(
        filename,
        samples_p,
        runs_p,
        output_tag,
        preprocess_min,
        max_len,
        flags,
        verbose=True
    )

    cd_hit_est_2d(
        formatted_fasta_fp,
        os.path.join(clusters_p, config["Filenames"]["cluster"]),
        formatted_fasta_base,
        float(config["CD-HIT"]["identity"]),
        logs_p,
        log_tag=output_tag,
        verbose=True
    )

    assignments = cluster_assign(
        f'format_{output_tag}.clstr',
        formatted_fasta_fp,
        runs_p
    )

    dict_cluster = json_load(os.path.join(clusters_p, config["Filenames"]["cluster_pkl"]))

    cluster_compile(
        assignments,
        dict_cluster,
        os.path.join(runs_p, f'format_{output_tag}.assign'),
        reports_p
    )

    cluster_report = cluster_miner(reports_p, output_tag, runs_p, runs_p, flags)

    threads = int(config['blast']['num_threads'])
    num_seqs = int(config['blast']['max_target_seqs'])

    if flags['Master']['C_BLAST']:
        blast_cluster = bclust(
            clusters_p,
            config["Filenames"]["cluster_metadata"],
            runs_p,
            blasts_p,
            config["Filenames"]["cluster"],
            runs_p,
            output_tag,
            flags,
            num_threads=threads,
            max_tar_seq=num_seqs,
            fasta=f"{output_tag}_to_blast.fasta"
        )

        if flags['Master']['L_BLAST']:
            blast_report = reblast(
                metadata_p,
                config["Filenames"]["metadata"],
                runs_p,
                blasts_p,
                config["Filenames"]["l_blast"],
                runs_p,
                output_tag,
                threads=threads,
                max_tar_seq=num_seqs,
                fasta=f"{output_tag}_to_reblast.fasta"
            )
            report_compiler(
                cluster_report,
                runs_p,
                output_tag,
                mappings,
                reports_p,
                flags,
                bclust_dict=blast_cluster,
                blast_dict=blast_report
            )
        else:
            report_compiler(
                cluster_report,
                runs_p,
                output_tag,
                mappings,
                reports_p,
                flags,
                bclust_dict=blast_cluster
            )
    else:
        report_compiler(
            cluster_report,
            runs_p,
            output_tag,
            mappings,
            reports_p,
            flags
        )

    if mode == 'contig' and single:
        id_report_fp = os.path.join(reports_p, f"{output_tag}_ID_Report.txt")
        best_segments, min_used = compute_best_segments_iterative(
            id_report_fp,
            formatted_fasta_fp,
            mappings,
            selection_start_min,
            step_bp=contig_single_min_step_bp,
            floor_bp=contig_single_min_floor_bp,
        )
        flags['Sample']['best_segments'] = best_segments
        flags['Sample']['best_segments_min_len'] = min_used

    redirector(
        f"{output_tag}_ID_Report.txt",
        flags,
        output_tag,
        mappings,
        runs_p,
        reports_p,
        force_flumut=flags['Master']['flumut'],
        force_genin=flags['Master']['genin'],
        force_getref=flags['Master']['getref'],
        mode=mode,
        single_sample=single
    )

    if flags['Master']['flumut']:
        conform_to_flumut(flags, runs_p, output_tag)
        if update_flumut.upper() == 'ON':
            update_flumut_db()
        run_flumut(reports_p, runs_p, output_tag)
        remap_flumut_report(reports_p, mappings, output_tag)
        mut_miner(os.path.join(reports_p, f'{output_tag}_markers.tsv'), muts_of_interest, flags)

    if flags['Master']['nextclade']:
        fasta_to_nextclade(flags, runs_p, output_tag)
        run_nextclade(f'{output_tag}_to_nextclade', flags, reports_p, runs_p, output_tag)

        if flags['Final Report']['Sequences for NextClade']['H1'] != []:
            remap_nextclade(reports_p, f"{output_tag}_H1_nextclade.tsv", mappings)
            flags['Sample']['clade'] = mine_clade(
                os.path.join(reports_p, f"{output_tag}_H1_nextclade.tsv"),
                flags,
                mappings
            )

        if flags['Final Report']['Sequences for NextClade']['H3'] != []:
            remap_nextclade(reports_p, f"{output_tag}_H3_nextclade.tsv", mappings)
            flags['Sample']['clade'] = mine_clade(
                os.path.join(reports_p, f"{output_tag}_H3_nextclade.tsv"),
                flags,
                mappings
            )

        if flags['Final Report']['Sequences for NextClade']['H5'] != []:
            remap_nextclade(reports_p, f"{output_tag}_H5_nextclade.tsv", mappings)
            flags['Sample']['clade'] = mine_clade(
                os.path.join(reports_p, f"{output_tag}_H5_nextclade.tsv"),
                flags,
                mappings
            )

    if flags['Master']['getref']:
        get_reference(
            flags['Final Report']['Get References'],
            references_p,
            config['Filenames']['ref_db'],
            output_tag,
            email=config['getref']['email'],
            update_local_db=True
        )

    if single or forced_genin:
        flags['Master']['genin'] = to_genin2(flags)
    else:
        flags['Master']['genin']= False

    if flags['Master']['genin']:
        run_genin2(runs_p, output_tag, os.path.join(reports_p, f'{output_tag}_genin.tsv'))
        remap_genin2(reports_p, f'{output_tag}_genin.tsv', mappings, flags)

    if mode == 'consensus':
        generate_final_report(
            os.path.join(reports_p, f"{output_tag}_ID_Report.txt"),
            flags,
            seg_lens,
            mappings,
            muts_loci_meaning,
            html_skeleton,
            os.path.join(reports_p, f'{output_tag}_final_report.html'),
            f'{output_tag}_final_report',
            runs_p=runs_p
            )
        with open(os.path.join(reports_p, f'{output_tag}_flags.json'), 'w') as f:
            json.dump(flags, f)

        print(f'Finished: {filename}')
        return {
            'flags': flags,
            'output_tag': output_tag
        }
def main(flagdict=flagdict):
    """Main entry point for the pipeline."""
    args = parser()
    cwd = os.getcwd()

    config_file = args.config
    mode = args.mode
    update_flumut = args.update_flumut_db
    off_apps = args.turn_off
    force_apps = args.force
    single = args.single_sample.lower() == 'on'

    config = configparser.ConfigParser()
    read_files = config.read(config_file)

    if not read_files:
        print(f"Could not read configuration file: {config_file}")
        sys.exit(1)

    required_sections = ['Paths', 'Functions', 'Sequence_Size', 'Filenames', 'blast', 'CD-HIT', 'getref']
    missing_sections = [sec for sec in required_sections if sec not in config]
    if missing_sections:
        print(f"Missing config sections: {', '.join(missing_sections)}")
        sys.exit(1)

    # Paths from config
    samples = config['Paths']['samples']
    runs = config['Paths']['runs']
    references = config['Paths']['references']
    reports = config['Paths']['reports']
    logs = config['Paths']['logs']
    blasts = config['Paths']['blast_database']
    clusters = config['Paths']['cluster_database']
    metadata = config['Paths']['metadata']

    rm_previous = config['Functions']['remove_previous'] if args.remove_previous.lower() == 'on' else False

    samples_p = os.path.abspath(os.path.join(cwd, samples))
    runs_p = os.path.abspath(os.path.join(cwd, runs))
    references_p = os.path.abspath(os.path.join(cwd, references))
    reports_p = os.path.abspath(os.path.join(cwd, reports))
    logs_p = os.path.abspath(os.path.join(cwd, logs))
    blasts_p = os.path.abspath(os.path.join(cwd, blasts))
    clusters_p = os.path.abspath(os.path.join(cwd, clusters))
    metadata_p = os.path.abspath(os.path.join(cwd, metadata))

    max_len = int(config["Sequence_Size"]["max"]) if args.max_length is None else args.max_length
    min_len = int(config["Sequence_Size"]["min"]) if args.min_length is None else args.min_length

    print('Configuration file:', config_file)
    print('Mode:', mode)
    print('Batch mode:', args.batch)
    print('Single sample mode:', single)
    print('Remove previous files:', rm_previous)
    print('Update Flumut DB:', update_flumut)
    print('Max sequence length:', max_len)
    print('Min sequence length:', min_len)
    print('Batch directory:', args.batch_dir if args.batch_dir else '.')

    # Check required paths
    for pth, label in [
        (samples_p, "Samples"),
        (runs_p, "Runs"),
        (references_p, "References"),
        (reports_p, "Reports"),
        (blasts_p, "Blast database"),
        (clusters_p, "Cluster database"),
        (metadata_p, "Metadata"),
    ]:
        if not os.path.exists(pth):
            print(f"{label} path does not exist: {pth}")
            sys.exit(1)

    if not os.path.exists(logs_p):
        os.makedirs(logs_p)

    # Optional cleanup
    if rm_previous:
        # Clean runs directory
        for file in glob.glob(os.path.join(runs_p, "*")):
            if os.path.isfile(file):
                os.remove(file)

        # Clean generated FASTA/helper files in samples directory root
        for pattern in [
            "*_to_blast.fasta",
            "*_to_reblast.fasta",
            "*_to_flumut.fasta",
            "*_to_nextclade_H1.fasta",
            "*_to_nextclade_H3.fasta",
            "*_to_nextclade_H5.fasta",
            "*_to_genin.fasta",
            "format_*"
        ]:
            for file in glob.glob(os.path.join(samples_p, pattern)):
                if os.path.isfile(file):
                    os.remove(file)

    # Batch mode
    if args.batch:
        try:
            fasta_files = find_fasta_files(samples_p, batch_dir=args.batch_dir, recursive=True)
        except Exception as e:
            print(f"Could not scan batch directory: {e}")
            sys.exit(1)

        if not fasta_files:
            batch_target = args.batch_dir if args.batch_dir else '.'
            print(f"No FASTA files found under samples/{batch_target}")
            sys.exit(1)

        print(f"Discovered {len(fasta_files)} FASTA files")

        failures = []
        results = []

        batch_label = args.batch_dir.strip().replace(os.sep, "_").replace("/", "_") if args.batch_dir else "batch_root"
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

        for i, fasta_file in enumerate(fasta_files, start=1):
            print(f"\n[{i}/{len(fasta_files)}] Processing {fasta_file}")
            try:
                result = run_pipeline_for_file(
                    filename=fasta_file,
                    flags_template=flagdict,
                    config=config,
                    samples_p=samples_p,
                    runs_p=runs_p,
                    references_p=references_p,
                    reports_p=reports_p,
                    logs_p=logs_p,
                    blasts_p=blasts_p,
                    clusters_p=clusters_p,
                    metadata_p=metadata_p,
                    mode=mode,
                    single=single,
                    update_flumut=update_flumut,
                    off_apps=off_apps,
                    force_apps=force_apps,
                    rm_previous=rm_previous,
                    max_len=max_len,
                    min_len=min_len,
                    batch=True,
                    batch_dir=batch_label,
                    timestamp=timestamp
                )

                results.append({
                    'input_file': fasta_file,
                    'file_tag': result['output_tag'],
                    'status': 'ok',
                    'error': ''
                })

            except Exception as e:
                print(f"ERROR processing {fasta_file}: {e}")
                failures.append((fasta_file, str(e)))

                sample_stem = Path(fasta_file).stem
                failed_tag = f"{batch_label}_{timestamp}_{sample_stem}"

                results.append({
                    'input_file': fasta_file,
                    'file_tag': failed_tag,
                    'status': 'failed',
                    'error': str(e)
                })

        batch_summary_fp = os.path.join(
            reports_p,
            f"batch_summary_{batch_label}_{timestamp}.tsv"
        )

        with open(batch_summary_fp, 'w') as out:
            out.write('input_file\tfile_tag\tstatus\terror\n')
            for row in results:
                out.write(f"{row['input_file']}\t{row['file_tag']}\t{row['status']}\t{row['error']}\n")

        print(f"\nBatch summary written to: {batch_summary_fp}")

        maybe_create_batch_artifacts(reports_p, batch_summary_fp)

        if failures:
            print("\nBatch completed with failures:")
            for fasta_file, err in failures:
                print(f" - {fasta_file}: {err}")
            sys.exit(1)

        print("\nBatch completed successfully")
        return

    # Single-file mode
    run_pipeline_for_file(
        filename=args.filename,
        flags_template=flagdict,
        config=config,
        samples_p=samples_p,
        runs_p=runs_p,
        references_p=references_p,
        reports_p=reports_p,
        logs_p=logs_p,
        blasts_p=blasts_p,
        clusters_p=clusters_p,
        metadata_p=metadata_p,
        mode=mode,
        single=single,
        update_flumut=update_flumut,
        off_apps=off_apps,
        force_apps=force_apps,
        rm_previous=rm_previous,
        max_len=max_len,
        min_len=min_len,
        batch=False,
        batch_dir='',
        timestamp=''
    )
if __name__=='__main__':
    main()
