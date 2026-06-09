#!/usr/bin/env python3
#### IMPORT DEPENDENCIES
import os
import sys
import argparse
import configparser
import subprocess
import json
import ast
import csv
import glob
import re
import hashlib
import shutil
import tempfile
import numpy as np
from typing import Any, Iterable, Optional
from collections import defaultdict
from flu_utils import seq_get,dict_to_fasta,headers_from_mult_fas,parse_clstr,json_load,\
    concat_fasta,\
    seq_filter_get,mine_genotype_H,mine_genotype_N,convert_to_prop, mine_single_HA,mine_single_NA
from metadata_utils import ClusterReportTable, ClusterMetadata, SequenceMetadata
from metadata_tree_utils import load_reference_dicts
from structures import flagdict, muts_loci_meaning , int_to_iupac, muts_interest,seg_lens, \
    iav_segments, contig_single_min_floor_bp, contig_single_min_step_bp
from final_report_utils import (
    html_skeleton,
    export_cluster_composition,
    export_id_report_rollup,
    generate_final_report,
    maybe_create_batch_artifacts,
    order_mutation_labels,
)
from copy import deepcopy
import pandas as pd
from gb_utils import fetch_genbank_list, append_genbank_from_list, get_gb, create_lookup,\
    write_gb,gb_to_fasta
from pathlib import Path
from datetime import datetime


_TOOL_VERSION_CACHE: dict[str, dict[str, str]] = {}

_BATCH_SUMMARY_FIELDS = [
    'input_file', 'file_tag', 'status', 'error', 'sample_dir', 'reports_dir_rel',
    'flumut_ran', 'nextclade_ran', 'genin2_ran', 'clade', 'genotype',
    'genin2_genotype', 'genin2_subgenotype',
    'genin2_pb2', 'genin2_pb1', 'genin2_pa', 'genin2_np', 'genin2_na', 'genin2_mp', 'genin2_ns',
    'mutations_of_interest',
]

_FLUMUT_BATCH_SEPARATOR = '__SEP__'


def _set_tool_status(flags: dict, tool: str, status: str, detail: str = '') -> None:
    """Store pipeline tool execution status in the shared flags object."""
    flags.setdefault('Tool Status', {})
    flags['Tool Status'][tool] = {
        'status': status,
        'detail': detail,
    }


def _write_flags_json(reports_p: str, output_tag: str, flags: dict) -> None:
    """Persist pipeline flags for post-run inspection."""
    with open(os.path.join(reports_p, f'{output_tag}_flags.json'), 'w') as f:
        json.dump(flags, f)


def _write_tool_versions_tsv(reports_p: str, output_tag: str, tool_versions: dict[str, str]) -> None:
    """Persist tool version metadata as a simple TSV."""
    output_fp = os.path.join(reports_p, f'{output_tag}_tool_versions.tsv')
    with open(output_fp, 'w', encoding='utf-8') as out:
        out.write('tool\tversion\n')
        for tool in sorted(tool_versions):
            version = str(tool_versions[tool]).replace('\t', ' ').replace('\n', ' ').strip()
            out.write(f'{tool}\t{version}\n')


def _version_command_for(tool_name: str) -> Optional[list[str]]:
    """Return the command used to query a CLI tool version."""
    commands = {
        'blastn': ['blastn', '-version'],
        'cd-hit-est-2d': ['cd-hit-est-2d', '-h'],
        'flumut': ['flumut', '-V'],
        'genin2': ['genin2', '--version'],
        'nextclade': ['nextclade', '--version'],
    }
    return commands.get(tool_name)


def _first_nonempty_line(text: str) -> str:
    """Return the first non-empty line from command output."""
    for line in text.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped
    return ''


def _parse_tool_versions(tool_name: str, output_text: str) -> dict[str, str]:
    """Normalize raw version command output into tool/version pairs."""
    output_text = str(output_text).strip()
    if not output_text:
        return {tool_name: 'unavailable'}

    if tool_name == 'flumut':
        parsed = {}
        for line in output_text.splitlines():
            stripped = line.strip()
            if stripped.startswith('FluMut:'):
                parsed['flumut'] = stripped.split(':', 1)[1].strip() or 'unavailable'
            elif stripped.startswith('FluMutDB:'):
                parsed['FluMutDB'] = stripped.split(':', 1)[1].strip() or 'unavailable'
        if parsed:
            parsed.setdefault('flumut', 'unavailable')
            parsed.setdefault('FluMutDB', 'unavailable')
            return parsed
        return {'flumut': _first_nonempty_line(output_text) or 'unavailable', 'FluMutDB': 'unavailable'}

    if tool_name == 'blastn':
        line = _first_nonempty_line(output_text)
        return {'blastn': line.split(':', 1)[1].strip() if ':' in line else (line or 'unavailable')}

    if tool_name == 'cd-hit-est-2d':
        match = re.search(r'CD-HIT version\s+([^\s]+)\s+\(built on\s+([^)]+)\)', output_text, flags=re.IGNORECASE)
        if match:
            return {'cd-hit-est-2d': f'{match.group(1)} (built on {match.group(2)})'}
        return {'cd-hit-est-2d': _first_nonempty_line(output_text) or 'unavailable'}

    return {tool_name: _first_nonempty_line(output_text) or 'unavailable'}


def _capture_tool_version(tool_name: str) -> dict[str, str]:
    """Capture and cache version metadata for a supported CLI tool."""
    if tool_name in _TOOL_VERSION_CACHE:
        return dict(_TOOL_VERSION_CACHE[tool_name])

    command = _version_command_for(tool_name)
    if command is None:
        versions = {tool_name: 'unavailable'}
    else:
        try:
            completed = subprocess.run(command, capture_output=True, text=True, check=False)
            output_text = '\n'.join(part for part in [completed.stdout, completed.stderr] if part).strip()
            versions = _parse_tool_versions(tool_name, output_text)
        except Exception:
            versions = {tool_name: 'unavailable'}

    _TOOL_VERSION_CACHE[tool_name] = dict(versions)
    return dict(versions)


def _record_tool_versions(command: list[str], tool_versions: Optional[dict[str, str]], enabled: bool = False) -> None:
    """Record CLI tool versions for the executable used by a command."""
    if not enabled or tool_versions is None or not command:
        return

    executable = os.path.basename(str(command[0]).strip())
    if not executable:
        return

    for key, value in _capture_tool_version(executable).items():
        tool_versions[key] = value


def _require_files(paths: Iterable[str], context: str) -> None:
    """Raise a clear error when expected pipeline artifacts are missing."""
    missing = [path for path in paths if not os.path.exists(path)]
    if missing:
        missing_str = ', '.join(missing)
        raise FileNotFoundError(f'{context}. Missing required file(s): {missing_str}')


def _run_command(command: list[str], tool_name: str, sample_tag: str = '', expected_outputs: Optional[Iterable[str]] = None,
                 stdout: Any = None, stderr: Any = None, tool_versions: Optional[dict[str, str]] = None,
                 capture_tool_versions: bool = False) -> None:
    """Run an external command and raise a contextual error on failure."""
    _record_tool_versions(command, tool_versions, enabled=capture_tool_versions)
    completed = subprocess.run(command, stdout=stdout, stderr=stderr)
    if completed.returncode != 0:
        sample_msg = f' for {sample_tag}' if sample_tag else ''
        raise RuntimeError(f'{tool_name} failed{sample_msg} with exit code {completed.returncode}')
    if expected_outputs:
        _require_files(expected_outputs, f'{tool_name} completed but did not produce all expected outputs')

### STEP FUNCTIONS
#### FASTA PREPROCESS
def seq_enum(counter: int) -> str:
    """
    Generate an internal sequence identifier.

    Parameters
    ----------
    counter : int
        Counter value to encode in the identifier.

    Returns
    -------
    str
        Sequence identifier in the ``>Seq-<counter>`` format.
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
) -> dict[str, str]:
    """
    Filter, rename, and export sequences for a pipeline run.

    Parameters
    ----------
    filename : str
        Relative FASTA path inside ``samples_path``.
    samples_path : str
        Directory containing the input FASTA file.
    runs_path : str
        Directory where the formatted FASTA file is written.
    file_tag : str
        Flat identifier used for generated filenames.
    min_seq_len : int
        Minimum accepted sequence length.
    max_seq_len : int
        Maximum accepted sequence length.
    flags : dict
        Pipeline state dictionary updated with rejected sequence headers.
    verbose : bool, default False
        Whether filtering decisions should be printed.

    Returns
    -------
    dict
        Mapping from internal sequence identifiers to original FASTA headers.
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

    if not fasta_dict:
        failure_msg = (
            f'No sequences remained after preprocessing for {filename} '
            f'after applying length thresholds {min_seq_len}-{max_seq_len}.'
        )
        flags.setdefault('Fasta Preprocess', {})
        flags['Fasta Preprocess']['Failure'] = failure_msg
        raise ValueError(failure_msg)

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
def cd_hit_est_2d(filename: str, cluster_reps: str, output: str, identity: float, logs_p: str, log_tag: str = None,
                  verbose: bool = False, tool_versions: Optional[dict[str, str]] = None,
                  capture_tool_versions: bool = False) -> None:
    """
    Run `cd-hit-est-2d` against cluster representatives.

    Parameters
    ----------
    filename : str
        Input FASTA file to classify.
    cluster_reps : str
        FASTA file of cluster representative sequences.
    output : str
        Output file stem for CD-HIT results.
    identity : float
        Sequence identity threshold.
    logs_p : str
        Directory where stdout and stderr logs are written when
        ``verbose=False``.
    log_tag : str, optional
        Filename stem to use for log files.
    verbose : bool, default False
        Whether to stream subprocess output directly to the terminal.

    Returns
    -------
    None
    """
    tag = log_tag if log_tag is not None else os.path.basename(filename)

    if not verbose:
        log_fp = os.path.join(logs_p, f'{tag}_cd_hit.log')
        err_fp = os.path.join(logs_p, f'{tag}_cd_hit.err')
        with open(log_fp, 'a') as out, open(err_fp, 'a') as err:
            _run_command(
                ['cd-hit-est-2d', '-i', cluster_reps, '-i2', filename, '-o', output, '-c', str(identity), '-g', '1'],
                'cd-hit-est-2d',
                sample_tag=tag,
                expected_outputs=[f'{output}.clstr', output],
                stdout=out,
                stderr=err,
                tool_versions=tool_versions,
                capture_tool_versions=capture_tool_versions,
            )
    else:
        _run_command(
            ['cd-hit-est-2d', '-i', cluster_reps, '-i2', filename, '-o', output, '-c', str(identity), '-g', '1'],
            'cd-hit-est-2d',
            sample_tag=tag,
            expected_outputs=[f'{output}.clstr', output],
            tool_versions=tool_versions,
            capture_tool_versions=capture_tool_versions,
        )

#### CLUSTERING AND CLUSTER REPORT
def cluster_assign(report:str, headers_file:str,runsdir:str) -> dict[str, str]:
    """
    Assign sample sequences to CD-HIT clusters.

    Parameters
    ----------
    report : str
        CD-HIT ``.clstr`` filename inside ``runsdir``.
    headers_file : str
        FASTA file whose headers should be assigned.
    runsdir : str
        Directory containing the cluster report and assignment output.

    Returns
    -------
    dict
        Mapping from normalized sample headers to cluster identifiers.

    Notes
    -----
    A ``.assign`` file is written to ``runsdir``.
    """
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

def cluster_compile(s_dict: dict[str, str], cl_dict: dict[str, list[dict[str, int] | str]], cl_assign:str, reports_path:str) -> None:
    """
    Combine cluster assignments with cluster metadata into a report.

    Parameters
    ----------
    s_dict : dict
        Sample-to-cluster assignment mapping.
    cl_dict : dict
        Cluster metadata lookup dictionary.
    cl_assign : str
        Path to the assignment file created by :func:`cluster_assign`.
    reports_path : str
        Directory where the compiled report is written.

    Returns
    -------
    None
    """
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
    
    #adding to the information sample:cluster,rep,identity the seg,gen, host, and collection date ranges
    for key in dict_sample:
        try:
            compiler[key]=[]
            compiler[key].append(sample_assign[key])
            compiler[key].append(dict_sample[key])
            if dict_sample[key]!='>Unassigned':
                representative, genotypes, segments, hosts, countries, collection_date = _unpack_cluster_payload(
                    dict_clust[dict_sample[key]]
                )
                compiler[key].extend([genotypes, segments, hosts, countries, collection_date, representative])
        except KeyError:
            print(f'Key: {key} not found')
            continue
    output_name=cl_assign.split('/')[-1]
    output_name=output_name.replace('.assign','')
    outpath=os.path.join(reports_path,f'{output_name.replace("format_","")}_clust_report.txt')
    
    #outputting report
    with open(f'{outpath}','w') as report:
        report.write('Sample_access\t%ID\tCluster\tCluster_rep\tGenotypes\tSegment\tHosts\tCountries\tCollection_date\n')
        for sample in compiler:
            if len(compiler[sample])>2:
                report.write(f'{sample}\t{compiler[sample][0]}\t\
                            {compiler[sample][1]}\t{compiler[sample][7]}\t\
                            {compiler[sample][2]}\t{compiler[sample][3]}\t{compiler[sample][4]}\t{compiler[sample][5]}\t{compiler[sample][6]}\n')
            elif len(compiler[sample])==0:
                continue  
            else:
                report.write(f'{sample}\t{compiler[sample][0]}\t{compiler[sample][1]}\n')
    print(f"Compiled report generated in {outpath}")

def cluster_miner(reports_p: str, file_tag: str, formatted_fasta_dir: str, runs_p: str, flags: dict) -> dict[str, list[str | dict[str, int | float]]]:
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
        safe_row += [""] * max(0, 8 - len(safe_row))

        to_report[key_str] = [
            str(safe_row[2]).replace(" ", "") if safe_row[2] is not None else "",
            str(safe_row[1]).replace("                        >", ">") if safe_row[1] is not None else "",
            str(safe_row[0]).replace("%", "") if safe_row[0] is not None else "",
            safe_row[4] if safe_row[4] is not None else "",
            safe_row[3] if safe_row[3] is not None else "",
            safe_row[5] if safe_row[5] is not None else "",
            safe_row[6] if safe_row[6] is not None else "",
            _normalize_collection_date(safe_row[7]),
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
def best_blast(runs_p:str,blast_report:str,thres_id:float=0.0) -> dict[str, list[str | float]]:
    """
    Select the best BLAST hit for each query sequence.

    Parameters
    ----------
    runs_p : str
        Directory containing the BLAST output file.
    blast_report : str
        BLAST tabular report filename.
    thres_id : float, default 0.0
        Minimum percent identity required to keep a hit.

    Returns
    -------
    dict
        Mapping from query sequence name to ``[representative, identity]``.
    """
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

def bclust(metadata_p:str,metadata_f:str,samples_p:str,blast_p:str,blast_db:str,runs_p:str,file_tag:str,flags:dict,
           num_threads:int=2,max_tar_seq:int=7,fasta:str='to_blast.fasta', tool_versions: Optional[dict[str, str]] = None,
           capture_tool_versions: bool = False) -> dict[str, list[str | float | dict[str, int | float]]]:
    """
    Run cluster-representative BLAST for unassigned sequences.

    Parameters
    ----------
    metadata_p : str
        Directory containing the cluster metadata table.
    metadata_f : str
        Cluster metadata filename.
    samples_p : str
        Directory containing the FASTA file to BLAST.
    blast_p : str
        Directory containing the BLAST database.
    blast_db : str
        BLAST database name.
    runs_p : str
        Directory where BLAST outputs are written.
    file_tag : str
        Flat sample identifier used to name output files.
    flags : dict
        Pipeline state dictionary updated with re-BLAST routing decisions.
    num_threads : int, default 2
        Number of BLAST threads.
    max_tar_seq : int, default 7
        Maximum number of target sequences to request.
    fasta : str, default 'to_blast.fasta'
        FASTA filename to analyze inside ``samples_p``.

    Returns
    -------
    dict
        Best-hit cluster annotations keyed by query sequence name.
    """
    access=headers_from_mult_fas([os.path.join(samples_p,fasta)],only_name=True)
    queries=[key for key in access]
    blast_output = os.path.join(runs_p, f"{file_tag}_bcrun.txt")
    _run_command(
        ['blastn', '-db', os.path.join(blast_p,blast_db), '-query', os.path.join(samples_p,fasta),
         '-out', blast_output, '-outfmt', '6', '-num_threads', str(num_threads), '-max_target_seqs', str(max_tar_seq)],
        'cluster BLAST',
        sample_tag=file_tag,
        expected_outputs=[blast_output],
        tool_versions=tool_versions,
        capture_tool_versions=capture_tool_versions,
    )
    b_tab=best_blast(runs_p,f"{file_tag}_bcrun.txt",thres_id=90.0)
    pd.set_option('display.max_colwidth', None)
    metadata=pd.read_table(os.path.join(metadata_p,metadata_f),index_col=False)
    metadata_columns = {str(col).strip().upper(): col for col in metadata.columns}
    cluster_col = metadata_columns.get('CLUSTER', 'Cluster')
    representative_col = metadata_columns.get('REPRESENTATIVE', 'Representative')
    segments_col = metadata_columns.get('SEGMENTS', 'Segments')
    genotypes_col = metadata_columns.get('GENOTYPES', 'Genotypes')
    hosts_col = metadata_columns.get('HOSTS', 'Hosts')
    countries_col = metadata_columns.get('COUNTRIES', 'Countries')
    collection_date_col = metadata_columns.get('COLLECTION_DATE') or metadata_columns.get('COLLECTION_DATE'.replace('_', ''))
    assigned=set(b_tab.keys())
    queries=set(queries)
    to_reblast=queries-assigned
    for key in b_tab:
        b_tab[key].append(metadata[cluster_col][metadata[representative_col]==b_tab[key][0]].to_string(index=False,min_rows=None,max_rows=None))
        b_tab[key].append(metadata[segments_col][metadata[representative_col]==b_tab[key][0]].to_string(index=False,min_rows=None,max_rows=None))
        b_tab[key][3]=ast.literal_eval(b_tab[key][3])
        b_tab[key].append(metadata[genotypes_col][metadata[representative_col]==b_tab[key][0]].to_string(index=False,min_rows=None,max_rows=None))
        b_tab[key][4]=ast.literal_eval(b_tab[key][4])
        b_tab[key].append(metadata[hosts_col][metadata[representative_col]==b_tab[key][0]].to_string(index=False,min_rows=None,max_rows=None))
        b_tab[key][5]=ast.literal_eval(b_tab[key][5])
        b_tab[key].append(metadata[countries_col][metadata[representative_col]==b_tab[key][0]].to_string(index=False,min_rows=None,max_rows=None))
        b_tab[key][6]=ast.literal_eval(b_tab[key][6])
        if collection_date_col is not None:
            date_value = metadata[collection_date_col][metadata[representative_col]==b_tab[key][0]].to_string(index=False,min_rows=None,max_rows=None)
        else:
            date_value = ''
        b_tab[key].append(_normalize_collection_date(date_value))

    if list(to_reblast) != []:
        for i in to_reblast:
            flags['BLAST']['Sequences unassigned against cluster representatives'].append(i)
        to_reblast_dict=seq_filter_get(os.path.join(samples_p,fasta),to_reblast)
        dict_to_fasta(to_reblast_dict,os.path.join(runs_p,f'{file_tag}_to_reblast'))
        flags['Master']['L_BLAST']=True
    else:
        flags['Master']['L_BLAST']=False
    return b_tab

def reblast(metadata_p:str,metadata_f:str,samples_p:str,blast_p:str,blast_db:str,runs_p:str,file_tag:str,threads:int=2,
            max_tar_seq:int=7,fasta:str='to_reblast.fasta', tool_versions: Optional[dict[str, str]] = None,
            capture_tool_versions: bool = False) -> dict[str, list[str | float]]:
    """
    Run local-database BLAST for sequences still unassigned after cluster BLAST.

    Parameters
    ----------
    metadata_p : str
        Directory containing the local metadata table.
    metadata_f : str
        Metadata filename.
    samples_p : str
        Directory containing the FASTA file to BLAST.
    blast_p : str
        Directory containing the BLAST database.
    blast_db : str
        BLAST database name.
    runs_p : str
        Directory where BLAST outputs are written.
    file_tag : str
        Flat sample identifier used to name output files.
    threads : int, default 2
        Number of BLAST threads.
    max_tar_seq : int, default 7
        Maximum number of target sequences to request.
    fasta : str, default 'to_reblast.fasta'
        FASTA filename to analyze inside ``samples_p``.

    Returns
    -------
    dict
        Best-hit local metadata annotations keyed by query sequence name.
    """
    blast_output = os.path.join(runs_p, f"{file_tag}_brun.txt")
    _run_command(
        ['blastn', '-db', os.path.join(blast_p,blast_db), '-query', os.path.join(samples_p,fasta),
         '-out', blast_output, '-outfmt', '6', '-num_threads', str(threads), '-max_target_seqs', str(max_tar_seq)],
        'local BLAST',
        sample_tag=file_tag,
        expected_outputs=[blast_output],
        tool_versions=tool_versions,
        capture_tool_versions=capture_tool_versions,
    )
    report=best_blast(runs_p,f"{file_tag}_brun.txt")
    metadata=pd.read_csv(os.path.join(metadata_p,metadata_f),sep=';',index_col=False)
    metadata_columns = {str(col).upper(): col for col in metadata.columns}
    required_columns = ['ACCESSION', 'SEGMENT', 'GENOTYPE', 'HOST', 'COUNTRY']
    missing_columns = [column for column in required_columns if column not in metadata_columns]
    if missing_columns:
        raise ValueError(
            f"Metadata file {os.path.join(metadata_p, metadata_f)} is missing required columns: {', '.join(missing_columns)}"
        )

    accession_col = metadata_columns['ACCESSION']
    segment_col = metadata_columns['SEGMENT']
    genotype_col = metadata_columns['GENOTYPE']
    host_col = metadata_columns['HOST']
    country_col = metadata_columns['COUNTRY']
    collection_date_col = metadata_columns.get('COLLECTION_DATE')
    for key in report:
        report[key].append(metadata[segment_col][metadata[accession_col]==report[key][0]].to_string(index=False))
        report[key].append(metadata[genotype_col][metadata[accession_col]==report[key][0]].to_string(index=False))
        report[key].append(metadata[host_col][metadata[accession_col]==report[key][0]].to_string(index=False))
        report[key].append(metadata[country_col][metadata[accession_col]==report[key][0]].to_string(index=False))
        if collection_date_col is not None:
            date_value = metadata[collection_date_col][metadata[accession_col]==report[key][0]].to_string(index=False)
            report[key].append(_normalize_collection_date(date_value))
        else:
            report[key].append('Unknown')
    return report

#### REPORT COMPILER
def report_compiler(clust_dict: dict[str, list[str | dict[str, int | float]]], samples_p: str, file_tag: str, mappings_dict: dict[str, str], reports_p: str, flags: dict, bclust_dict: Optional[dict[str, list[str | float | dict[str, int | float]]]] = {}, blast_dict: Optional[dict[str, list[str | float]]] = {}) -> None:
    """
    Merge cluster and BLAST assignments into the final ID report.

    Parameters
    ----------
    clust_dict : dict
        Cluster-derived report entries.
    samples_p : str
        Directory containing the formatted sample FASTA.
    file_tag : str
        Flat sample identifier used to name output files.
    mappings_dict : dict
        Mapping from internal FASTA identifiers to original sequence headers.
    reports_p : str
        Directory where the final report is written.
    flags : dict
        Pipeline state dictionary used during report generation.
    bclust_dict : dict, default {}
        Cluster-BLAST report entries.
    blast_dict : dict, default {}
        Local BLAST report entries.

    Returns
    -------
    None
    """
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
                to_report[key].append(_normalize_collection_date(bclust_dict[key][7]))
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
                to_report[key].append(_normalize_collection_date(blast_dict[key][6]))
                to_report[key].append('L-BLAST')
    assigned=set(to_report.keys())
    seqs=set(seqs)
    to_remote=seqs-assigned
    mappings=mappings_dict
    if list(to_remote)!=[]:
        for seq in to_remote:
            to_report[seq]=['Unassigned','NA','NA','NA','NA','NA','NA','Unknown','NA']
    with open(os.path.join(reports_p, f'{file_tag}_ID_Report.txt'), 'w') as report:
        report.write('Sample_name\tSegment\tCluster\tRepresentative\t%ID\tAssigned_by\tGenotype\tHost\tCountry\tCollection_date\n'.upper())
        for key in to_report:
            mapped=mappings[f'>{key}']
            #print(to_report[key])
            report.write(f'{mapped}\t{to_report[key][3]}\t\
                         {to_report[key][1]}\t{to_report[key][0]}\t{to_report[key][2]}\t\
                            {to_report[key][8]}\t{to_report[key][4]}\t{to_report[key][5]}\t{to_report[key][6]}\t{to_report[key][7]}\n')
    if len(list(to_remote))>0:
        flags['BLAST']['Sequences unassigned against local database'].extend(to_remote)
    #print(f'Report generated in {os.path.join(reports_p,filename.replace(".fasta",""))}_ID_Report.txt')

#### REPORT PARSERS
def _parse_segment_cell(x: Any) -> Optional[int]:
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


def _parse_pident_cell(x: Any) -> Optional[float]:
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


def compute_best_segments_iterative(id_report_fp: str, formatted_fasta_fp: str, mappings: dict[str, str], start_min_len: int,
                                   step_bp: int = contig_single_min_step_bp,
                                   floor_bp: int = contig_single_min_floor_bp) -> tuple:
    """
    Select a single best contig per segment without overwriting selections.

    Parameters
    ----------
    id_report_fp : str
        Path to the ID report TSV.
    formatted_fasta_fp : str
        Path to the formatted FASTA file.
    mappings : dict
        Mapping from internal FASTA identifiers to original sequence headers.
    start_min_len : int
        Initial minimum contig length threshold.
    step_bp : int, default ``contig_single_min_step_bp``
        Step size used when relaxing the minimum contig length.
    floor_bp : int, default ``contig_single_min_floor_bp``
        Lowest minimum contig length allowed during iterative relaxation.

    Returns
    -------
    tuple
        ``(best_segments_internal, min_len_used)`` where the first item maps
        canonical segment names to selected internal sequence identifiers.

    Notes
    -----
    Candidates are ranked by contig length and then by percent identity. Once a
    segment is filled it is not replaced in later iterations.
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

    def _len_for_internal(sid: Any) -> int:
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
def safe_literal_eval(x: Any) -> Any:
    """
    Safely parse a scalar or container represented as a string.

    Parameters
    ----------
    x : Any
        Value to parse. Empty values and values that cannot be parsed are
        returned unchanged.

    Returns
    -------
    Any
        Parsed Python literal when parsing succeeds, otherwise the original
        value.
    """
    if pd.isna(x):
        return x
    s = str(x).strip()
    if not s:
        return x
    try:
        return ast.literal_eval(s)
    except Exception:
        return x


def _normalize_collection_date(value: Any) -> str:
    """Normalize missing or blank collection dates for report output."""
    if pd.isna(value):
        return 'Unknown'

    text = str(value).strip()
    return text if text else 'Unknown'


def _unpack_cluster_payload(payload: list[Any]) -> tuple[str, Any, Any, Any, Any, str]:
    """Support both legacy and Collection_date-aware cluster payload shapes."""
    if len(payload) >= 6:
        return (
            str(payload[5]).strip(),
            payload[0],
            payload[1],
            payload[2],
            payload[3],
            _normalize_collection_date(payload[4]),
        )

    return (
        str(payload[4]).strip() if len(payload) >= 5 else '',
        payload[0] if len(payload) >= 1 else {},
        payload[1] if len(payload) >= 2 else {},
        payload[2] if len(payload) >= 3 else {},
        payload[3] if len(payload) >= 4 else {},
        'Unknown',
    )
def redirector(report: str, flags: dict, file_tag: str, mappings: dict[str, str], runs_p: str, reports_p: str,
               force_flumut: bool, force_genin: bool, force_getref: bool, mode: str, single_sample: bool = True) -> None:
    """
    Populate downstream analysis targets from an ID report.

    Parameters
    ----------
    report : str
        Name of the ID report file inside ``reports_p``.
    flags : dict
        Mutable pipeline state dictionary updated in place.
    file_tag : str
        Flat sample identifier used in generated filenames.
    mappings : dict
        Mapping from internal FASTA identifiers to original sequence headers.
    runs_p : str
        Path to the run workspace.
    reports_p : str
        Path to the reports directory.
    force_flumut : bool
        Whether FluMut should be forced even when routing rules would skip it.
    force_genin : bool
        Whether GenIn2 should be forced.
    force_getref : bool
        Whether reference retrieval should be forced.
    mode : {"contig", "consensus"}
        Pipeline execution mode.
    single_sample : bool, default True
        Whether the current run is treated as a single-sample analysis.

    Returns
    -------
    None

    Notes
    -----
    This function mutates ``flags`` in place, fills the ``Sample`` and
    ``Final Report`` sections, and may rewrite the report file when contig mode
    is restricted to the selected best segments.
    """
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
    if 'COLLECTION_DATE' in report_df.columns:
        report_df['COLLECTION_DATE'] = report_df['COLLECTION_DATE'].apply(_normalize_collection_date)
    report_df=report_df.dropna(subset=['SAMPLE_NAME', 'SEGMENT'])
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
def update_flumut_db(tool_versions: Optional[dict[str, str]] = None, capture_tool_versions: bool = False) -> None:
    """Update the local FluMut database using the external CLI."""
    _run_command(
        ['flumut', '--update'],
        'FluMut DB update',
        tool_versions=tool_versions,
        capture_tool_versions=capture_tool_versions,
    )

def conform_to_flumut(
    flagdict: dict,
    sample_path: str,
    file_tag: str,
    mappings: Optional[dict[str, str]] = None,
    separator: str = '_',
    use_original_headers: bool = False,
    output_suffix: str = 'to_flumut',
) -> str:
    """
    Write a FluMut-ready FASTA containing only selected sequences.

    Parameters
    ----------
    flagdict : dict
        Pipeline state dictionary containing per-segment assignments and the
        list of sequences selected for FluMut.
    sample_path : str
        Directory containing ``format_{file_tag}.fasta``.
    file_tag : str
        Flat sample identifier used to name the output FASTA.

    Returns
    -------
    str
        Output FASTA path.

    Notes
    -----
    Output headers are suffixed with the segment name so they match FluMut's
    expected naming convention.
    """
    seq_seg = {}
    for key in flagdict['Sample']:
        if key in ('PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS'):
            for i in flagdict['Sample'][key]:
                if i in flagdict['Final Report']["Sequences for FluMut"]:
                    seq_seg[i] = key
    flt_fasta = seq_get(os.path.join(sample_path, f'format_{file_tag}.fasta'))
    outdict = {}
    used_headers = set()

    if mappings is None:
        mappings = {}

    for key in flt_fasta:
        if key in seq_seg:
            segment = seq_seg[key]
            header_base = key
            if use_original_headers:
                header_base = str(mappings.get(key, key)).strip()
                if header_base.startswith('>'):
                    header_base = header_base[1:]
                if not header_base:
                    header_base = key.lstrip('>')

            output_header = f'>{header_base}{separator}{segment}'
            if output_header in used_headers:
                duplicate_index = 2
                while True:
                    candidate = f'>{header_base}__DUP{duplicate_index}{separator}{segment}'
                    if candidate not in used_headers:
                        output_header = candidate
                        break
                    duplicate_index += 1

            used_headers.add(output_header)
            outdict[output_header] = flt_fasta[key]

    output_base = os.path.join(sample_path, f'{file_tag}_{output_suffix}')
    dict_to_fasta(outdict, output_base)
    return f'{output_base}.fasta'

def conform_to_genin(flagdict: dict, mappings: dict[str, str], sample_path: str, file_tag: str) -> None:
    """
    Write a GenIn2-ready FASTA preserving a shared sample prefix.

    Parameters
    ----------
    flagdict : dict
        Pipeline state dictionary containing per-segment assignments and the
        list of sequences selected for GenIn2.
    mappings : dict
        Mapping from internal FASTA identifiers to original sequence headers.
    sample_path : str
        Directory containing ``format_{file_tag}.fasta``.
    file_tag : str
        Flat sample identifier used to name the output FASTA.

    Returns
    -------
    None

    Notes
    -----
    GenIn2 groups records by the header stem before the final ``_SEGMENT``
    suffix, so this export uses the sample file stem across all segments.
    """
    seq_seg = {}
    for key in flagdict['Sample']:
        if key in ('PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS'):
            for header in flagdict['Sample'][key]:
                if header in flagdict['Final Report']['Sequences for GenIn']:
                    seq_seg[header] = key

    formatted_fasta = seq_get(os.path.join(sample_path, f'format_{file_tag}.fasta'))
    outdict = {}

    for selected_header, segment in seq_seg.items():
        internal_header = selected_header if selected_header in formatted_fasta else None
        if internal_header is None and selected_header in mappings:
            internal_header = selected_header

        if internal_header is None or internal_header not in formatted_fasta:
            continue

        outdict[f'>{file_tag}_{segment}'] = formatted_fasta[internal_header]

    dict_to_fasta(outdict, os.path.join(sample_path, f'{file_tag}_to_genin'))
        
def run_flumut(
    reports_p: str,
    samples_p: str,
    file_tag: str,
    regex: str = "(.+)\_(.+)",
    excel_output: Optional[str] = None,
    write_tsv_outputs: bool = True,
    input_fasta: Optional[str] = None,
    tool_versions: Optional[dict[str, str]] = None,
    capture_tool_versions: bool = False,
) -> None:
    """
    Run FluMut on a prepared FASTA.

    Parameters
    ----------
    reports_p : str
        Directory where FluMut report files are written.
    samples_p : str
        Directory containing the FluMut-ready FASTA file.
    file_tag : str
        Flat sample identifier used to build input and output filenames.
    regex : str, default ``"(.+)\_(.+)"``
        Regular expression passed to FluMut for parsing FASTA headers.
    excel_output : str, optional
        Excel workbook path passed to FluMut via ``-x``.
    write_tsv_outputs : bool, default True
        Whether to request the standard TSV outputs.
    input_fasta : str, optional
        Explicit FASTA path to run instead of ``{file_tag}_to_flumut.fasta``.

    Returns
    -------
    None
    """
    print("Running flumut")
    markers_fp = os.path.join(reports_p, f"{file_tag}_markers.tsv")
    mutations_fp = os.path.join(reports_p, f"{file_tag}_mutations.tsv")
    literature_fp = os.path.join(reports_p, f"{file_tag}_literature.tsv")
    fasta_fp = input_fasta if input_fasta is not None else os.path.join(samples_p, f"{file_tag}_to_flumut.fasta")
    command = ['flumut']
    expected_outputs = []

    if write_tsv_outputs:
        command.extend(['-m', markers_fp, '-M', mutations_fp, '-l', literature_fp])
        expected_outputs.extend([markers_fp, mutations_fp])

    if excel_output:
        command.extend(['-x', excel_output])
        expected_outputs.append(excel_output)

    command.extend([fasta_fp, '-n', regex])

    _run_command(
        command,
        'FluMut',
        sample_tag=file_tag,
        expected_outputs=expected_outputs,
        tool_versions=tool_versions,
        capture_tool_versions=capture_tool_versions,
    )
def remap_flumut_report(reports_dir: str, mappings: dict, filename: str) -> None:
    """
    Replace internal sequence identifiers in FluMut reports.

    Parameters
    ----------
    reports_dir : str
        Directory containing the FluMut output tables.
    mappings : dict
        Mapping from internal FASTA identifiers to original sequence headers.
    filename : str
        File stem used by the FluMut output files.

    Returns
    -------
    None

    Notes
    -----
    The markers and mutations TSV files are rewritten in place.
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
def mut_miner(dataframe: str, muts_of_interest: dict[str, Iterable[str]], flagsdict: dict) -> None:
    """
    Extract mutations of interest from a FluMut report.

    Parameters
    ----------
    dataframe : str
        Path to the FluMut markers report TSV file.
    muts_of_interest : dict
        Mapping of segment names to the mutations that should be tracked.
    flagsdict : dict
        Pipeline state dictionary updated with per-segment mutation hits.

    Returns
    -------
    None

    Notes
    -----
    Matching mutations are appended to ``flagsdict['Sample']["<SEG>_muts"]``.
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
def get_reference(listref: Iterable[str], references_p: str, reference_output_p: str, db_p: str, file_tag: str, email: Optional[str] = None, update_local_db: bool = False) -> None:
    """
    Retrieve reference GenBank records and export them for a sample.

    Parameters
    ----------
    listref : list
        Reference accession identifiers to retrieve.
    references_p : str
        Directory containing the shared local reference database.
    reference_output_p : str
        Directory where the sample-specific reference files are written.
    db_p : str
        Filename of the local GenBank database inside ``references_p``.
    file_tag : str
        Flat sample identifier used to name the output files.
    email : str, optional
        Email address used by NCBI Entrez when remote retrieval is needed.
    update_local_db : bool, default False
        Whether remotely fetched records should also be appended to the local
        reference database.

    Returns
    -------
    None
    """
    references=set(listref)
    found=[]
    missing=[]
    reference_db_fp = os.path.join(references_p, db_p)
    output_p = reference_output_p

    #checking if local_db has entries
    records=create_lookup(reference_db_fp)
    for ref in references:
        if ref in records:
            found.append(ref)
        else:
            missing.append(ref)
    #getting missing records
    remote=fetch_genbank_list(missing,email)
    #getting local records
    local=get_gb(reference_db_fp, found)
    #creating output files
    output=remote+local
    #writing output files (.gb and .fasta)
    outname = f'{file_tag}_references'
    write_gb(output,os.path.join(output_p,f'{outname}.gb'))
    gb_to_fasta(os.path.join(output_p,f'{outname}.gb'),os.path.join(output_p,f'{outname}.fasta'))

    if update_local_db:
        append_genbank_from_list(reference_db_fp,remote)
#### NEXTCLADE PATH
def fasta_to_nextclade(flagsdict: dict, samples_p: str, file_tag: str) -> None:
    """
    Split selected HA sequences into Nextclade input FASTA files.

    Parameters
    ----------
    flagsdict : dict
        Pipeline state dictionary containing the sequences routed to Nextclade.
    samples_p : str
        Directory containing ``format_{file_tag}.fasta``.
    file_tag : str
        Flat sample identifier used to name the generated FASTA files.

    Returns
    -------
    None
    """
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
def run_nextclade(fasta: str,flagsdict: dict,reports_p: str,samples_p: str,file_tag: str,
                  tool_versions: Optional[dict[str, str]] = None, capture_tool_versions: bool = False) -> None:
    """
    Run Nextclade for each supported HA genotype present in the sample.

    Parameters
    ----------
    fasta : str
        Prefix of the genotype-specific FASTA files to analyze.
    flagsdict : dict
        Pipeline state dictionary containing the sequences selected for
        Nextclade.
    reports_p : str
        Directory where Nextclade TSV reports are written.
    samples_p : str
        Directory containing the genotype-specific FASTA files.
    file_tag : str
        Flat sample identifier used to name output reports.

    Returns
    -------
    None
    """
    dict_builds_broad={'H1':'flu_h1n1pdm_ha_broad','H3':'flu_h3n2_ha_broad','H5':'community/moncla-lab/iav-h5/ha/all-clades'}
    for key in flagsdict['Final Report']['Sequences for NextClade']:
        if flagsdict['Final Report']['Sequences for NextClade'][key] != []:
            fasta_gen=key
            nextclade_fp = os.path.join(reports_p, f"{file_tag}_{fasta_gen}_nextclade.tsv")
            _run_command(
                ['nextclade', 'run', '-d', dict_builds_broad[fasta_gen], '--output-tsv', nextclade_fp, os.path.join(samples_p,f"{fasta}_{fasta_gen}.fasta")],
                f'Nextclade {fasta_gen}',
                sample_tag=file_tag,
                expected_outputs=[nextclade_fp],
                tool_versions=tool_versions,
                capture_tool_versions=capture_tool_versions,
            )
def remap_nextclade(reports_p: str,nextclade_report: str,mappings_dict: dict[str, str]) -> None:
    """
    Replace internal sequence identifiers in a Nextclade report.

    Parameters
    ----------
    reports_p : str
        Directory containing the Nextclade TSV file.
    nextclade_report : str
        Filename of the Nextclade TSV report.
    mappings_dict : dict
        Mapping from internal FASTA identifiers to original sequence headers.

    Returns
    -------
    None
    """
    nextclade_rep=pd.read_table(os.path.join(reports_p,nextclade_report),index_col=False)
    nextclade_rep['seqName']=nextclade_rep['seqName'].apply(lambda x: mappings_dict[f'>{x}'])
    nextclade_rep.to_csv(os.path.join(reports_p,nextclade_report), sep='\t', index=False)

#### GENIN PATH

def mine_clade(clade_file: str, flagsdict: dict, mappings_dict: dict[str, str]) -> str:
    """
    Extract the HA clade assignment from a Nextclade report.

    Parameters
    ----------
    clade_file : str
        Path to the Nextclade TSV report.
    flagsdict : dict
        Pipeline state dictionary containing HA genotype and selected sequence
        information.
    mappings_dict : dict
        Mapping from internal FASTA identifiers to original sequence headers.

    Returns
    -------
    str
        Clade assignment, ``"Unassigned"``, ``"Unassigned by error"``, or
        ``"Non-applicable"`` depending on the report contents and routing
        context.
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

def to_genin2(flagsdict: dict) -> bool:
    """
    Decide whether the sample should be sent to GenIn2.

    Parameters
    ----------
    flagsdict : dict
        Pipeline state dictionary containing the inferred HA genotype and
        clade.

    Returns
    -------
    bool
        True when the sample is H5 and belongs to clade ``2.3.4.4b``,
        otherwise False.
    """
    if flagsdict["Sample"]["H_gen"] == "H5":
        if flagsdict["Sample"]["clade"] == "2.3.4.4b":
            return True
    return False

def run_genin2(samples_p: str, file_tag: str, genin_report: str,
               tool_versions: Optional[dict[str, str]] = None, capture_tool_versions: bool = False) -> None:
    """
    Run GenIn2 on the FluMut-formatted sample FASTA.

    Parameters
    ----------
    samples_p : str
        Directory containing ``{file_tag}_to_genin.fasta``.
    file_tag : str
        Flat sample identifier used to build the input FASTA name.
    genin_report : str
        Output directory passed to GenIn2.

    Returns
    -------
    None
    """
    _run_command(
        ['genin2','-o', genin_report, os.path.join(samples_p,f"{file_tag}_to_genin.fasta")],
        'GenIn2',
        sample_tag=file_tag,
        expected_outputs=[genin_report],
        tool_versions=tool_versions,
        capture_tool_versions=capture_tool_versions,
    )

def remap_genin2(reports_p: str, genin_report: str, mappings: dict[str, str], flagsdict: dict) -> dict:
    """
    Rewrite GenIn2 sample names and store parsed constellation results.

    Parameters
    ----------
    reports_p : str
        Directory containing the GenIn2 report.
    genin_report : str
        Filename of the GenIn2 TSV report.
    mappings : dict
        Mapping from internal FASTA identifiers to original sequence headers.
    flagsdict : dict
        Pipeline state dictionary updated with parsed GenIn2 genotype data.

    Returns
    -------
    dict
        The updated ``flagsdict``.

    Notes
    -----
    The input report is rewritten in place after sample names are remapped.
    """
    genin_path = os.path.join(reports_p, genin_report)
    genin_df = pd.read_table(genin_path, index_col=False)

    segment_cols = ['PB2', 'PB1', 'PA', 'NP', 'NA', 'MP', 'NS']
    segment_to_original = {}
    for segment in ('PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS'):
        originals = []
        for internal_header in flagsdict.get('Sample', {}).get(segment, []):
            original_header = mappings.get(internal_header)
            if original_header:
                originals.append(str(original_header).strip())
        if originals:
            segment_to_original[segment] = originals

    original_stems = set()
    for segment, originals in segment_to_original.items():
        suffix = f'_{segment}'
        for original in originals:
            header_text = original.lstrip('>')
            if header_text.endswith(suffix):
                original_stems.add(header_text[:-len(suffix)])

    shared_original_stem = list(original_stems)[0] if len(original_stems) == 1 else ''

    def _remap_genin_sample_name(sample_name: Any, row: pd.Series) -> str:
        name = str(sample_name).strip()
        valid = []
        for segment in segment_cols:
            value = row.get(segment)
            if pd.notna(value) and str(value).strip() not in ('', '?'):
                valid.append(segment)

        if len(valid) == 1 and valid[0] in segment_to_original and len(segment_to_original[valid[0]]) == 1:
            return segment_to_original[valid[0]][0]

        if shared_original_stem:
            return shared_original_stem

        return name

    genin_df['Sample Name'] = [
        _remap_genin_sample_name(sample_name, row)
        for sample_name, (_, row) in zip(genin_df['Sample Name'].tolist(), genin_df.iterrows())
    ]

    genin_df.to_csv(genin_path, sep='\t', index=False)

    result = []
    skipped = []

    for _, row in genin_df.iterrows():
        sample_name = str(row.get('Sample Name', '')).strip()
        genotype = str(row.get('Genotype', '')).strip()

        sub_genotype = row.get('Sub-genotype', '')
        if pd.isna(sub_genotype) or str(sub_genotype).strip() == '':
            sub_genotype = 'N/A'
        else:
            sub_genotype = str(sub_genotype).strip()

        valid = row[segment_cols][row[segment_cols].notna() & (row[segment_cols].astype(str).str.strip() != '?')]

        if len(valid) == 0:
            skipped.append({'sample': sample_name, 'reason': 'no_valid_segment'})
            continue

        row_result = {
            'Sample Name': sample_name,
            'Genotype': genotype,
            'Sub-genotype': sub_genotype,
            'Notes': '' if pd.isna(row.get('Notes', '')) else str(row.get('Notes', '')).strip(),
        }
        for segment in segment_cols:
            value = row.get(segment, '')
            row_result[segment] = '' if pd.isna(value) else str(value).strip()

        if len(valid) > 1:
            row_result['Output Format'] = 'constellation'
        else:
            row_result['Output Format'] = 'expanded'

        result.append(row_result)

    flagsdict['Sample']['Genin_genotypes'] = result
    flagsdict['Sample']['Genin_skipped_rows'] = skipped
    return flagsdict


#### MAIN
def parser() -> argparse.Namespace:
    """
    Parse command-line arguments for the main pipeline entry point.

    Returns
    -------
    argparse.Namespace
        Parsed command-line options for the AFluID pipeline.
    """
    parser = argparse.ArgumentParser(
        prog='AFluID',
        description='AFluID: Automated Influenza Identification Pipeline'
    )
    parser.add_argument('-c', '--config', type=str, help='Configuration file path', default=os.path.join(os.path.dirname(__file__), 'config.ini'), required=False)
    parser.add_argument('-f', '--filename', type=str, help='FASTA file name', required=False)
    parser.add_argument('-m', '--mode', type=str, help='Mode of operation', choices=('contig', 'consensus'), required=True)
    parser.add_argument('-b', '--batch', action='store_true', help='Run pipeline for all FASTA files in batch directory')
    parser.add_argument('-bd', '--batch_dir', type=str, default='', help='Subdirectory inside samples/ to use for batch mode')
    parser.add_argument('-bn', '--batch_name', type=str, help='Optional name for the batch output directory', required=False)
    parser.add_argument('-ff', '--force', help='Force run of aditional tools', nargs='*', choices=('flumut', 'genin', 'getref'))
    parser.add_argument('-fdb', '--update_flumut_db', type=str, help='turn off auto-update for flumut db', choices=('on', 'off'), default='on', required=False)
    parser.add_argument('-ss', '--single_sample', type=str, help='Single sample mode', default='on', choices=('on', 'off'), required=False)
    parser.add_argument('-o', '--output_name', type=str, help='Optional output stem for direct non-batch single-file runs', required=False)
    parser.add_argument('-off', '--turn_off', help='Turn Off additional analysis tools', nargs='*', choices=('flumut', 'genin', 'nextclade', 'getref'))
    parser.add_argument('-rm', '--remove_previous', type=str, help='Remove previous files', default='on', choices=('on', 'off'), required=False)
    parser.add_argument('-Ml', '--max_length', type=int, help='Maximum sequence length', required=False)
    parser.add_argument('-ml', '--min_length', type=int, help='Minimum sequence length', required=False)
    parser.add_argument('--outdir', type=str, help='Optional root directory for per-sample report/reference outputs', required=False)
    parser.add_argument('--replace', action='store_true', help='Replace existing output directories before running')
    parser.add_argument('--skip_cdhit', action='store_true', help='Bypass cd-hit-est-2d and send all sequences directly to cluster BLAST')
    parser.add_argument('-tv', '--tool-versions', action='store_true', help='Write an optional tool/version TSV from always-captured CLI version metadata')
    parser.add_argument('-bf', '--batch_fasta', action='store_true', help='Demultiplex one multifasta into a temporary batch of per-sample FASTAs')
    parser.add_argument('--force-bf', action='store_true', help='If --batch_fasta input is not demultipliable, run it as a multi-sample FASTA instead of failing')
    parser.add_argument('--bf-regex', type=str, required=False, help='Mandatory regex for --batch_fasta that captures only the sample identifier')

    args = parser.parse_args()

    if not args.batch and not args.filename:
        parser.error("--filename is required unless --batch is used")

    if args.batch_name and not (args.batch or args.batch_fasta):
        parser.error("--batch_name can only be used with --batch or --batch_fasta")

    if args.batch_name is not None and not str(args.batch_name).strip():
        parser.error("--batch_name cannot be empty")

    if args.batch and args.batch_fasta:
        parser.error("--batch and --batch_fasta cannot be used together")

    if args.batch_fasta and not args.filename:
        parser.error("--batch_fasta requires --filename")

    if args.batch_fasta and not args.bf_regex:
        parser.error("--bf-regex is required when --batch_fasta is used")

    if args.bf_regex and not args.batch_fasta:
        parser.error("--bf-regex can only be used with --batch_fasta")

    if args.force_bf and not args.batch_fasta:
        parser.error("--force-bf can only be used with --batch_fasta")

    if args.output_name is not None:
        output_name = str(args.output_name).strip()
        if not output_name:
            parser.error("--output_name cannot be empty")
        if output_name in {'.', '..'} or '/' in output_name or '\\' in output_name:
            parser.error("--output_name must be a single name without path separators")
        if re.fullmatch(r'[A-Za-z0-9._-]+', output_name) is None:
            parser.error("--output_name may only contain letters, numbers, dots, underscores, and hyphens")
        if args.batch or args.batch_fasta:
            parser.error("--output_name can only be used in direct non-batch single-file runs")

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


def resolve_output_root(cwd: str, outdir: Optional[str]) -> Optional[str]:
    """Resolve an optional output root from cwd unless it is absolute."""
    if not outdir:
        return None

    path = Path(outdir).expanduser()
    if not path.is_absolute():
        path = Path(cwd) / path

    return str(path.resolve())


def resolve_config_path(raw_path: str, project_root: Path) -> str:
    """Resolve config-managed paths relative to the project root unless absolute."""
    path = Path(raw_path).expanduser()
    if path.is_absolute():
        return str(path.resolve())
    return str((project_root / path).resolve())


def resolve_single_sample_output_paths(
    filename: str,
    default_reports_p: str,
    default_references_p: str,
    custom_output_root: Optional[str] = None,
    output_name: Optional[str] = None,
) -> dict[str, str]:
    """Resolve per-sample report/reference output directories."""
    sample_folder = output_name if output_name else make_output_stem(filename, batch=False)
    sample_root = ''

    if custom_output_root:
        sample_root = os.path.join(custom_output_root, sample_folder)
        sample_reports_p = os.path.join(sample_root, 'reports')
        sample_reference_outputs_p = os.path.join(sample_root, 'references')
    else:
        sample_reports_p = os.path.join(default_reports_p, sample_folder)
        sample_reference_outputs_p = os.path.join(default_references_p, sample_folder)

    return {
        'sample_folder': sample_folder,
        'sample_root': sample_root if custom_output_root else '',
        'sample_reports_p': sample_reports_p,
        'sample_reference_outputs_p': sample_reference_outputs_p,
    }


def resolve_batch_name(batch_dir: str, timestamp: str, user_batch_name: Optional[str] = None) -> str:
    """Resolve the batch output directory name."""
    if user_batch_name is None:
        return make_run_prefix(batch_dir=batch_dir, timestamp=timestamp)

    batch_name = user_batch_name.strip()
    if not batch_name:
        raise ValueError('Batch name cannot be empty')

    if batch_name in {'.', '..'} or '/' in batch_name or '\\' in batch_name:
        raise ValueError('Batch name must be a single directory name without path separators')

    if re.fullmatch(r'[A-Za-z0-9._-]+', batch_name) is None:
        raise ValueError('Batch name may only contain letters, numbers, dots, underscores, and hyphens')

    return batch_name


def resolve_batch_output_paths(
    filename: str,
    default_reports_p: str,
    default_references_p: str,
    batch_name: str,
    custom_output_root: Optional[str] = None,
) -> dict[str, str]:
    """Resolve nested batch output directories for one sample."""
    sample_folder = make_batch_sample_tag(filename)

    if custom_output_root:
        batch_root = os.path.join(custom_output_root, batch_name)
        sample_root = os.path.join(batch_root, sample_folder)
        sample_reports_p = os.path.join(sample_root, 'reports')
        sample_reference_outputs_p = os.path.join(sample_root, 'references')
        reports_dir_rel = os.path.join(sample_folder, 'reports')
    else:
        batch_reports_root = os.path.join(default_reports_p, batch_name)
        batch_references_root = os.path.join(default_references_p, batch_name)
        sample_reports_p = os.path.join(batch_reports_root, sample_folder)
        sample_reference_outputs_p = os.path.join(batch_references_root, sample_folder)
        batch_root = batch_reports_root
        sample_root = sample_reports_p
        reports_dir_rel = sample_folder

    return {
        'batch_root': batch_root,
        'sample_folder': sample_folder,
        'sample_root': sample_root,
        'sample_reports_p': sample_reports_p,
        'sample_reference_outputs_p': sample_reference_outputs_p,
        'reports_dir_rel': reports_dir_rel,
    }


def run_batch_pipeline(
    fasta_files: list[str],
    batch_dir_label: str,
    search_root_label: str,
    batch_name: str,
    config: configparser.ConfigParser,
    batch_samples_p: str,
    runs_p: str,
    references_p: str,
    reports_p: str,
    logs_p: str,
    blasts_p: str,
    clusters_p: str,
    metadata_p: str,
    custom_output_root: Optional[str],
    replace_outputs: bool,
    mode: str,
    single: bool,
    update_flumut: str,
    off_apps: list[str] | None,
    force_apps: list[str] | None,
    rm_previous: bool,
    max_len: int,
    min_len: int,
    geo: str,
    taxa:str ,
    skip_cdhit: bool,
    write_tool_versions_tsv: bool,
) -> int:
    """Run the existing batch pipeline flow for a prepared FASTA list."""
    if not fasta_files:
        print(f"No FASTA files found under {search_root_label}")
        return 1

    print(f"Discovered {len(fasta_files)} FASTA files")

    failures = []
    results = []
    successful_flumut_fastas = []
    successful_batch_flumut_fastas = []
    batch_tool_versions = {}
    batch_flumut_error = ''

    print('Batch output name:', batch_name)

    if custom_output_root:
        batch_prepare_targets = [os.path.join(custom_output_root, batch_name)]
        batch_artifact_root = os.path.join(custom_output_root, batch_name)
    else:
        batch_prepare_targets = [
            os.path.join(reports_p, batch_name),
            os.path.join(references_p, batch_name),
        ]
        batch_artifact_root = os.path.join(reports_p, batch_name)

    try:
        prepare_output_targets(batch_prepare_targets, replace_outputs)
    except FileExistsError as exc:
        print(str(exc))
        return 1

    for i, fasta_file in enumerate(fasta_files, start=1):
        print(f"\n[{i}/{len(fasta_files)}] Processing {fasta_file}")
        batch_output_paths = resolve_batch_output_paths(
            fasta_file,
            reports_p,
            references_p,
            batch_name,
            custom_output_root,
        )
        try:
            result = run_pipeline_for_file(
                filename=fasta_file,
                flags_template=flagdict,
                config=config,
                samples_p=batch_samples_p,
                runs_p=runs_p,
                references_p=references_p,
                reports_p=batch_output_paths['sample_reports_p'],
                logs_p=logs_p,
                blasts_p=blasts_p,
                clusters_p=clusters_p,
                metadata_p=metadata_p,
                reference_output_p=batch_output_paths['sample_reference_outputs_p'],
                mode=mode,
                single=single,
                update_flumut=update_flumut,
                off_apps=off_apps,
                force_apps=force_apps,
                rm_previous=rm_previous,
                max_len=max_len,
                min_len=min_len,
                geo=geo,
                taxa=taxa,
                skip_cdhit=skip_cdhit,
                batch=True,
                batch_dir=batch_name,
                timestamp='',
                capture_tool_versions=True,
                write_tool_versions_tsv=False,
                batch_tool_versions=batch_tool_versions,
            )

            if not isinstance(result, dict) or 'output_tag' not in result:
                raise RuntimeError(
                    f"Pipeline completed without returning output metadata for {fasta_file}"
                )

            results.append({
                'input_file': fasta_file,
                'file_tag': result['output_tag'],
                'status': 'ok',
                'error': '',
                'sample_dir': batch_output_paths['sample_folder'],
                'reports_dir_rel': batch_output_paths['reports_dir_rel'],
                **extract_batch_summary_fields(result.get('flags', {})),
            })

            flumut_input_fp = str(result.get('flumut_input_fp', '')).strip()
            if flumut_input_fp:
                successful_flumut_fastas.append(flumut_input_fp)
            batch_flumut_input_fp = str(result.get('batch_flumut_input_fp', '')).strip()
            if batch_flumut_input_fp:
                successful_batch_flumut_fastas.append(batch_flumut_input_fp)

        except Exception as e:
            print(f"ERROR processing {fasta_file}: {e}")
            failures.append((fasta_file, str(e)))

            failed_tag = make_output_stem(
                fasta_file,
                batch=True,
                batch_dir=batch_name,
                timestamp='',
            )

            results.append({
                'input_file': fasta_file,
                'file_tag': failed_tag,
                'status': 'failed',
                'error': str(e),
                'sample_dir': batch_output_paths['sample_folder'],
                'reports_dir_rel': batch_output_paths['reports_dir_rel'],
                **{field: '' for field in _BATCH_SUMMARY_FIELDS if field not in ('input_file', 'file_tag', 'status', 'error', 'sample_dir', 'reports_dir_rel', 'flumut_ran', 'nextclade_ran', 'genin2_ran')},
                'flumut_ran': 'no',
                'nextclade_ran': 'no',
                'genin2_ran': 'no',
            })

    batch_summary_fp = os.path.join(batch_artifact_root, f"batch_summary_{batch_name}.tsv")

    with open(batch_summary_fp, 'w', encoding='utf-8', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=_BATCH_SUMMARY_FIELDS, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for row in results:
            normalized = {field: row.get(field, '') for field in _BATCH_SUMMARY_FIELDS}
            writer.writerow(normalized)

    print(f"\nBatch summary written to: {batch_summary_fp}")

    if successful_batch_flumut_fastas:
        batch_fasta_base = os.path.join(batch_artifact_root, f'{batch_name}_to_flumut_batch_tmp')
        batch_fasta_fp = f'{batch_fasta_base}.fasta'
        batch_excel_fp = os.path.join(batch_artifact_root, f'{batch_name}_flumut.xlsx')

        try:
            concat_fasta(successful_batch_flumut_fastas, batch_fasta_base)
            run_flumut(
                reports_p=batch_artifact_root,
                samples_p=runs_p,
                file_tag=batch_name,
                regex=f'(.+){_FLUMUT_BATCH_SEPARATOR}(.+)',
                excel_output=batch_excel_fp,
                write_tsv_outputs=False,
                input_fasta=batch_fasta_fp,
                tool_versions=batch_tool_versions,
                capture_tool_versions=True,
            )
        except Exception as exc:
            batch_flumut_error = str(exc)
            print(f'Warning: batch-wide FluMut workbook generation failed: {exc}')
        finally:
            if os.path.exists(batch_fasta_fp):
                os.remove(batch_fasta_fp)
    else:
        print('Skipping batch-wide FluMut workbook: no successful _to_flumut.fasta inputs were produced.')

    if write_tool_versions_tsv and batch_tool_versions:
        _write_tool_versions_tsv(batch_artifact_root, batch_name, batch_tool_versions)

    batch_cluster_composition_name = f'{batch_name}_cluster_composition.tsv'
    successful_id_report_fps = [
        os.path.join(
            batch_artifact_root,
            str(row.get('reports_dir_rel', '')).strip(),
            f"{str(row.get('file_tag', '')).strip()}_ID_Report.txt",
        )
        for row in results
        if str(row.get('status', '')).strip().lower() == 'ok' and str(row.get('file_tag', '')).strip()
    ]
    export_cluster_composition(
        id_report_fps=successful_id_report_fps,
        metadata_csv_fp=os.path.join(metadata_p, config['Filenames']['metadata']),
        cluster_clstr_fp=os.path.join(clusters_p, config['Filenames']['cluster_clstr']),
        output_fp=os.path.join(batch_artifact_root, batch_cluster_composition_name),
    )

    maybe_create_batch_artifacts(
        batch_artifact_root,
        batch_summary_fp,
        extra_files=[batch_cluster_composition_name],
    )

    if failures or batch_flumut_error:
        print("\nBatch completed with failures:")
        for fasta_file, err in failures:
            print(f" - {fasta_file}: {err}")
        if batch_flumut_error:
            print(f" - batch-wide FluMut workbook: {batch_flumut_error}")
        return 1

    print("\nBatch completed successfully")
    return 0


def prepare_output_targets(targets: list[str], replace: bool) -> None:
    """Create output directories, deleting existing ones only when requested."""
    seen = set()
    for target in targets:
        if not target or target in seen:
            continue
        seen.add(target)

        if os.path.exists(target):
            if not replace:
                raise FileExistsError(f'Output directory already exists: {target}. Use --replace to overwrite it.')

            if os.path.isdir(target):
                shutil.rmtree(target)
            else:
                raise FileExistsError(f'Output target exists and is not a directory: {target}')

        os.makedirs(target, exist_ok=True)


def create_run_workspace(base_runs_p: str) -> str:
    """Create a unique per-invocation intermediate workspace under runs/."""
    os.makedirs(base_runs_p, exist_ok=True)
    return tempfile.mkdtemp(prefix='session_', dir=base_runs_p)


def _tool_ran_yes_no(flags: dict, tool_name: str) -> str:
    """Return ``yes`` when a tool completed, otherwise ``no``."""
    status = str(flags.get('Tool Status', {}).get(tool_name, {}).get('status', '')).strip().lower()
    return 'yes' if status == 'completed' else 'no'


def _normalize_genin_summary(genin_rows: Any) -> dict[str, str]:
    """Collapse parsed GenIn2 rows into one batch-summary record."""
    summary = {
        'genin2_genotype': '',
        'genin2_subgenotype': '',
        'genin2_pb2': '',
        'genin2_pb1': '',
        'genin2_pa': '',
        'genin2_np': '',
        'genin2_na': '',
        'genin2_mp': '',
        'genin2_ns': '',
    }
    if not isinstance(genin_rows, list) or not genin_rows:
        return summary

    segment_map = {
        'PB2': 'genin2_pb2',
        'PB1': 'genin2_pb1',
        'PA': 'genin2_pa',
        'NP': 'genin2_np',
        'NA': 'genin2_na',
        'MP': 'genin2_mp',
        'NS': 'genin2_ns',
    }

    def _row_score(row: dict[str, Any]) -> int:
        return sum(1 for segment in segment_map if str(row.get(segment, '')).strip() not in ('', '?'))

    best_row = None
    best_score = -1
    for row in genin_rows:
        if not isinstance(row, dict):
            continue
        score = _row_score(row)
        if score > best_score:
            best_row = row
            best_score = score

    if best_row is None:
        return summary

    summary['genin2_genotype'] = str(best_row.get('Genotype', '')).strip()
    summary['genin2_subgenotype'] = str(best_row.get('Sub-genotype', '')).strip()

    for row in genin_rows:
        if not isinstance(row, dict):
            continue
        if not summary['genin2_genotype']:
            summary['genin2_genotype'] = str(row.get('Genotype', '')).strip()
        if not summary['genin2_subgenotype']:
            summary['genin2_subgenotype'] = str(row.get('Sub-genotype', '')).strip()
        for segment, field in segment_map.items():
            value = str(row.get(segment, '')).strip()
            if value in ('', '?') or summary[field]:
                continue
            summary[field] = value

    return summary


def _format_mutations_of_interest(flags: dict) -> str:
    """Aggregate all tracked mutations of interest into one semicolon-separated field."""
    mutations = []
    for segment in ('HA', 'NA', 'PA', 'PB1', 'PB2', 'MP', 'NP', 'NS'):
        for mut in flags.get('Sample', {}).get(f'{segment}_muts', []):
            mut = str(mut).strip()
            if not mut:
                continue
            label = muts_loci_meaning.get(mut, (f'{segment}:{mut}', ''))[0]
            mutations.append(label)
    return ';'.join(order_mutation_labels(mutations))


def extract_batch_summary_fields(flags: dict) -> dict[str, str]:
    """Build the enriched batch-summary fields from one sample's flags."""
    sample = flags.get('Sample', {})
    summary = {
        'flumut_ran': _tool_ran_yes_no(flags, 'flumut'),
        'nextclade_ran': _tool_ran_yes_no(flags, 'nextclade'),
        'genin2_ran': _tool_ran_yes_no(flags, 'genin'),
        'clade': str(sample.get('clade', '')).strip(),
        'genotype': str(sample.get('Genotype', '')).strip(),
        'mutations_of_interest': _format_mutations_of_interest(flags),
    }
    summary.update(_normalize_genin_summary(sample.get('Genin_genotypes', [])))
    return summary


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
        H5N1_IRMA_20260420_161518_PT_A_H5N1_SAMPLE_01_IRMA_a1b2c3d4e5
    """
    sample_stem = make_sample_stem(filename)

    if batch:
        run_prefix = make_run_prefix(batch_dir=batch_dir, timestamp=timestamp)
        sample_tag = make_batch_sample_tag(filename)
        return f"{run_prefix}_{sample_tag}"

    return sample_stem


def make_batch_sample_tag(filename: str, suffix_len: int = 10) -> str:
    """
    Build a readable per-sample batch tag.

    Files directly under ``samples/`` keep their basename. Files inside
    subdirectories append a short deterministic suffix derived from the
    relative parent path.
    """
    sample_stem = make_sample_stem(filename)
    parent = str(Path(filename).parent)

    if parent in ('', '.'):
        return sample_stem

    suffix = hashlib.sha1(parent.encode('utf-8')).hexdigest()[:suffix_len]
    return f"{sample_stem}_{suffix}"


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


def find_fasta_files(samples_p: str, batch_dir: str = '', recursive: bool = True) -> list[str]:
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


def resolve_single_sample_file(samples_p: str, filename: str) -> str:
    """
    Resolve a single requested FASTA file inside ``samples_p``.

    The input may be a relative path under ``samples/`` or a basename that is
    unique somewhere under that tree. The resolved return value is always a
    relative path from ``samples_p``.
    """
    requested = Path(filename)
    samples_root = Path(samples_p).resolve()
    direct = (samples_root / requested).resolve()

    if direct.exists() and direct.is_file():
        return str(direct.relative_to(samples_root))

    matches = sorted(
        p.relative_to(samples_root)
        for p in samples_root.rglob(requested.name)
        if p.is_file()
    )

    if not matches:
        raise FileNotFoundError(
            f"No sample file found under samples/ matching '{filename}'"
        )

    if len(matches) > 1:
        match_list = ', '.join(str(match) for match in matches)
        raise ValueError(
            f"Ambiguous sample filename '{filename}'. Matches: {match_list}"
        )

    return str(matches[0])


def _compile_batch_fasta_regex(pattern: str) -> re.Pattern[str]:
    """Compile and validate the mandatory demultiplex regex."""
    try:
        compiled = re.compile(pattern)
    except re.error as exc:
        raise ValueError(f"Invalid --bf-regex pattern: {exc}") from exc

    if 'sample' not in compiled.groupindex and compiled.groups < 1:
        raise ValueError("--bf-regex must define a named group 'sample' or at least one capture group")

    return compiled


def _extract_batch_fasta_sample_id(header: str, pattern: re.Pattern[str]) -> Optional[str]:
    """Extract the demultiplex sample identifier from a FASTA header."""
    match = pattern.search(header)
    if match is None:
        return None

    if 'sample' in pattern.groupindex:
        sample_id = match.group('sample')
    else:
        sample_id = match.group(1)

    sample_id = '' if sample_id is None else str(sample_id).strip()
    return sample_id or None


def _safe_batch_fasta_filename(sample_id: str, used_names: set[str]) -> str:
    """Convert a sample identifier into a filesystem-safe FASTA basename."""
    sanitized = sample_id.replace(os.sep, '_').replace('/', '_').strip()
    sanitized = re.sub(r'[^A-Za-z0-9._-]+', '_', sanitized)
    sanitized = sanitized.strip('._') or 'sample'

    candidate = sanitized
    if candidate in used_names:
        suffix = hashlib.sha1(sample_id.encode('utf-8')).hexdigest()[:10]
        candidate = f'{sanitized}_{suffix}'

    used_names.add(candidate)
    return candidate


def demultiplex_batch_fasta(input_fasta: str, output_dir: str, regex_pattern: str) -> dict[str, Any]:
    """Split one multifasta into per-sample FASTAs using a sample-ID regex."""
    compiled = _compile_batch_fasta_regex(regex_pattern)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    grouped_records: defaultdict[str, list[tuple[str, str]]] = defaultdict(list)
    unmatched_headers: list[str] = []
    current_header: Optional[str] = None
    current_sequence_lines: list[str] = []

    def save_current_record() -> None:
        nonlocal current_header, current_sequence_lines
        if current_header is None:
            return

        sample_id = _extract_batch_fasta_sample_id(current_header, compiled)
        if sample_id is None:
            unmatched_headers.append(current_header)
        else:
            sequence = ''.join(current_sequence_lines).replace(' ', '').strip()
            grouped_records[sample_id].append((current_header, sequence))

    with open(input_fasta, 'r', encoding='utf-8') as infile:
        for raw_line in infile:
            line = raw_line.rstrip('\n')
            if not line:
                continue

            if line.startswith('>'):
                save_current_record()
                current_header = line[1:].strip()
                current_sequence_lines = []
            else:
                current_sequence_lines.append(line.strip())

        save_current_record()

    demultipliable = len(grouped_records) > 0 and len(unmatched_headers) == 0
    written_files: list[str] = []
    sample_to_file: dict[str, str] = {}

    if demultipliable:
        used_names: set[str] = set()
        for sample_id in sorted(grouped_records):
            safe_name = _safe_batch_fasta_filename(sample_id, used_names)
            output_fp = output_path / f'{safe_name}.fasta'
            with open(output_fp, 'w', encoding='utf-8') as out:
                for header, sequence in grouped_records[sample_id]:
                    out.write(f'>{header}\n')
                    for i in range(0, len(sequence), 80):
                        out.write(f'{sequence[i:i + 80]}\n')
            written_files.append(str(output_fp))
            sample_to_file[sample_id] = str(output_fp)

    return {
        'demultipliable': demultipliable,
        'sample_ids': sorted(grouped_records.keys()),
        'unmatched_headers': unmatched_headers,
        'generated_fastas': written_files,
        'sample_to_file': sample_to_file,
    }

def run_pipeline_for_file(
    filename: str,
    flags_template: dict,
    config: configparser.ConfigParser,
    samples_p: str,
    runs_p: str,
    references_p: str,
    reports_p: str,
    logs_p: str,
    blasts_p: str,
    clusters_p: str,
    metadata_p: str,
    reference_output_p: Optional[str],
    mode: str,
    single: bool,
    update_flumut: str,
    off_apps: list[str] | None,
    force_apps: list[str] | None,
    rm_previous: bool,
    max_len: int,
    min_len: int,
    geo: str,
    taxa: str,
    skip_cdhit: bool = False,
    batch: bool = False,
    batch_dir: str = '',
    timestamp: str = '',
    output_name: Optional[str] = None,
    capture_tool_versions: bool = True,
    write_tool_versions_tsv: bool = False,
    batch_tool_versions: Optional[dict[str, str]] = None
) -> None:
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
    if batch_tool_versions is not None:
        tool_versions = batch_tool_versions
        flags['Tool Versions'] = dict(tool_versions)
    else:
        flags.setdefault('Tool Versions', {})
        tool_versions = flags['Tool Versions']


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
    print('Skip CD-HIT:', skip_cdhit)

    output_tag = make_output_stem(
        filename,
        batch=batch,
        batch_dir=batch_dir,
        timestamp=timestamp,
    )
    if not batch and output_name:
        output_tag = output_name

    active_reports_p = reports_p
    active_reference_output_p = reference_output_p if reference_output_p is not None else references_p

    os.makedirs(active_reports_p, exist_ok=True)
    os.makedirs(active_reference_output_p, exist_ok=True)

    print('Output tag:', output_tag)

    formatted_fasta_base = os.path.join(runs_p, f'format_{output_tag}')
    formatted_fasta_fp = f'{formatted_fasta_base}.fasta'

    preprocess_min = min_len
    if mode == 'contig' and single:
        preprocess_min = contig_single_min_floor_bp

    try:
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
    except Exception as exc:
        _set_tool_status(flags, 'preprocess', 'failed', str(exc))
        flags['Tool Versions'] = dict(tool_versions)
        _write_flags_json(active_reports_p, output_tag, flags)
        raise

    if not mappings:
        raise ValueError(
            f'No sequences remained after preprocessing for {filename}. '
            'Adjust the sequence length thresholds or inspect the input FASTA.'
        )

    _require_files([formatted_fasta_fp], 'FASTA preprocessing did not produce the formatted sample FASTA')
    _set_tool_status(flags, 'preprocess', 'completed')

    threads = int(config['blast']['num_threads'])
    num_seqs = int(config['blast']['max_target_seqs'])

    if skip_cdhit:
        cluster_report = {}
        flags['Master']['C_BLAST'] = True
        _set_tool_status(flags, 'cd_hit', 'skipped', 'Pipeline was run with --skip_cdhit')
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
            fasta=f"format_{output_tag}.fasta",
            tool_versions=tool_versions,
            capture_tool_versions=capture_tool_versions,
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
                fasta=f"{output_tag}_to_reblast.fasta",
                tool_versions=tool_versions,
                capture_tool_versions=capture_tool_versions,
            )
            report_compiler(
                cluster_report,
                runs_p,
                output_tag,
                mappings,
                active_reports_p,
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
                active_reports_p,
                flags,
                bclust_dict=blast_cluster
            )
    else:
        cd_hit_est_2d(
            formatted_fasta_fp,
            os.path.join(clusters_p, config["Filenames"]["cluster"]),
            formatted_fasta_base,
            float(config["CD-HIT"]["identity"]),
            logs_p,
            log_tag=output_tag,
            verbose=True,
            tool_versions=tool_versions,
            capture_tool_versions=capture_tool_versions,
        )
        _set_tool_status(flags, 'cd_hit', 'completed')

        assignments = cluster_assign(
            f'format_{output_tag}.clstr',
            formatted_fasta_fp,
            runs_p
        )
        _require_files(
            [os.path.join(runs_p, f'format_{output_tag}.clstr'), os.path.join(runs_p, f'format_{output_tag}.assign')],
            'CD-HIT clustering did not produce the expected assignment artifacts'
        )

        dict_cluster = json_load(os.path.join(clusters_p, config["Filenames"]["cluster_pkl"]))

        cluster_compile(
            assignments,
            dict_cluster,
            os.path.join(runs_p, f'format_{output_tag}.assign'),
            active_reports_p
        )

        cluster_report = cluster_miner(active_reports_p, output_tag, runs_p, runs_p, flags)

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
                fasta=f"{output_tag}_to_blast.fasta",
                tool_versions=tool_versions,
                capture_tool_versions=capture_tool_versions,
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
                    fasta=f"{output_tag}_to_reblast.fasta",
                    tool_versions=tool_versions,
                    capture_tool_versions=capture_tool_versions,
                )
                report_compiler(
                    cluster_report,
                    runs_p,
                    output_tag,
                    mappings,
                    active_reports_p,
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
                    active_reports_p,
                    flags,
                    bclust_dict=blast_cluster
                )
        else:
            report_compiler(
                cluster_report,
                runs_p,
                output_tag,
                mappings,
                active_reports_p,
                flags
            )

    id_report_fp = os.path.join(active_reports_p, f"{output_tag}_ID_Report.txt")
    _require_files([id_report_fp], 'Pipeline identification path did not produce the ID report')

    if mode == 'contig' and single:
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
        active_reports_p,
        force_flumut=flags['Master']['flumut'],
        force_genin=flags['Master']['genin'],
        force_getref=flags['Master']['getref'],
        mode=mode,
        single_sample=single
    )

    geo_dict, taxa_dict = load_reference_dicts(
        data_dir=metadata_p,
        geo_json_relpath=geo,
        tax_json_relpath=taxa,
    )
    export_id_report_rollup(
        id_report_fp=id_report_fp,
        output_fp=os.path.join(active_reports_p, f'{output_tag}_ID_Report_rollup.tsv'),
        geo_dict=geo_dict,
        taxa_dict=taxa_dict,
    )
    if not batch:
        export_cluster_composition(
            id_report_fps=[id_report_fp],
            metadata_csv_fp=os.path.join(metadata_p, config['Filenames']['metadata']),
            cluster_clstr_fp=os.path.join(clusters_p, config['Filenames']['cluster_clstr']),
            output_fp=os.path.join(active_reports_p, f'{output_tag}_cluster_composition.tsv'),
        )

    _set_tool_status(flags, 'flumut', 'skipped', 'FluMut routing not enabled for this sample')
    flumut_input_fp = os.path.join(runs_p, f'{output_tag}_to_flumut.fasta')
    batch_flumut_input_fp = ''
    if flags['Master']['flumut']:
        try:
            markers_fp = os.path.join(active_reports_p, f'{output_tag}_markers.tsv')
            mutations_fp = os.path.join(active_reports_p, f'{output_tag}_mutations.tsv')
            excel_output_fp = None if batch else os.path.join(active_reports_p, f'{output_tag}_flumut.xlsx')
            conform_to_flumut(flags, runs_p, output_tag)
            _require_files([flumut_input_fp], 'FluMut preparation did not produce the expected input FASTA')
            if batch:
                batch_flumut_input_fp = conform_to_flumut(
                    flags,
                    runs_p,
                    output_tag,
                    mappings=mappings,
                    separator=_FLUMUT_BATCH_SEPARATOR,
                    use_original_headers=True,
                    output_suffix='to_flumut_batch',
                )
                _require_files([batch_flumut_input_fp], 'Batch FluMut preparation did not produce the expected input FASTA')
            if update_flumut.upper() == 'ON':
                update_flumut_db(tool_versions=tool_versions, capture_tool_versions=capture_tool_versions)
            run_flumut(
                active_reports_p,
                runs_p,
                output_tag,
                excel_output=excel_output_fp,
                tool_versions=tool_versions,
                capture_tool_versions=capture_tool_versions,
            )
            _require_files([markers_fp, mutations_fp], 'FluMut did not produce the expected report files')
            remap_flumut_report(active_reports_p, mappings, output_tag)
            mut_miner(markers_fp, muts_of_interest, flags)
            _set_tool_status(flags, 'flumut', 'completed')
        except Exception as exc:
            flags['Master']['flumut'] = False
            flags['Final Report']['Sequences for FluMut'] = []
            flags['Final Report']['Sequences for GenIn'] = []
            _set_tool_status(flags, 'flumut', 'failed', str(exc))
            print(f'Warning: skipping FluMut for {output_tag}: {exc}')

    _set_tool_status(flags, 'nextclade', 'skipped', 'Nextclade routing not enabled for this sample')
    if flags['Master']['nextclade']:
        try:
            fasta_to_nextclade(flags, runs_p, output_tag)
            run_nextclade(
                f'{output_tag}_to_nextclade',
                flags,
                active_reports_p,
                runs_p,
                output_tag,
                tool_versions=tool_versions,
                capture_tool_versions=capture_tool_versions,
            )

            if flags['Final Report']['Sequences for NextClade']['H1'] != []:
                remap_nextclade(active_reports_p, f"{output_tag}_H1_nextclade.tsv", mappings)
                flags['Sample']['clade'] = mine_clade(
                    os.path.join(active_reports_p, f"{output_tag}_H1_nextclade.tsv"),
                    flags,
                    mappings
                )

            if flags['Final Report']['Sequences for NextClade']['H3'] != []:
                remap_nextclade(active_reports_p, f"{output_tag}_H3_nextclade.tsv", mappings)
                flags['Sample']['clade'] = mine_clade(
                    os.path.join(active_reports_p, f"{output_tag}_H3_nextclade.tsv"),
                    flags,
                    mappings
                )

            if flags['Final Report']['Sequences for NextClade']['H5'] != []:
                remap_nextclade(active_reports_p, f"{output_tag}_H5_nextclade.tsv", mappings)
                flags['Sample']['clade'] = mine_clade(
                    os.path.join(active_reports_p, f"{output_tag}_H5_nextclade.tsv"),
                    flags,
                    mappings
                )
            _set_tool_status(flags, 'nextclade', 'completed')
        except Exception as exc:
            flags['Master']['nextclade'] = False
            flags['Sample']['clade'] = ''
            flags['Final Report']['Sequences for NextClade'] = {'H1': [], 'H3': [], 'H5': []}
            _set_tool_status(flags, 'nextclade', 'failed', str(exc))
            print(f'Warning: skipping Nextclade for {output_tag}: {exc}')

    if flags['Master']['getref']:
        get_reference(
            flags['Final Report']['Get References'],
            references_p,
            active_reference_output_p,
            config['Filenames']['ref_db'],
            output_tag,
            email=config['getref']['email'],
            update_local_db=True
        )

    if forced_genin:
        flags['Master']['genin'] = True
    elif single:
        flags['Master']['genin'] = to_genin2(flags)
    else:
        flags['Master']['genin']= False

    if flags['Master']['genin']:
        flags['Final Report']['Sequences for GenIn'] = list(
            dict.fromkeys(flags['Final Report']['Sequences for FluMut'])
        )
    else:
        flags['Final Report']['Sequences for GenIn'] = []

    _set_tool_status(flags, 'genin', 'skipped', 'GenIn2 routing not enabled for this sample')
    if flags['Master']['genin']:
        try:
            genin_fp = os.path.join(active_reports_p, f'{output_tag}_genin.tsv')
            genin_input_fp = os.path.join(runs_p, f'{output_tag}_to_genin.fasta')
            conform_to_genin(flags, mappings, runs_p, output_tag)
            _require_files([genin_input_fp], 'GenIn2 preparation did not produce the expected input FASTA')
            run_genin2(
                runs_p,
                output_tag,
                genin_fp,
                tool_versions=tool_versions,
                capture_tool_versions=capture_tool_versions,
            )
            _require_files([genin_fp], 'GenIn2 did not produce the expected report file')
            remap_genin2(active_reports_p, f'{output_tag}_genin.tsv', mappings, flags)
            _set_tool_status(flags, 'genin', 'completed')
        except Exception as exc:
            flags['Master']['genin'] = False
            flags['Final Report']['Sequences for GenIn'] = []
            flags['Sample']['Genin_genotypes'] = []
            _set_tool_status(flags, 'genin', 'failed', str(exc))
            print(f'Warning: skipping GenIn2 for {output_tag}: {exc}')

    flags['Tool Versions'] = dict(tool_versions)
    _write_flags_json(active_reports_p, output_tag, flags)

    if mode == 'consensus':
        _require_files([id_report_fp], 'Final report generation requires the ID report')
        try:
            generate_final_report(
                        id_report_fp,
                        flags,
                        seg_lens,
                        mappings,
                        muts_loci_meaning,
                        html_skeleton,
                        os.path.join(active_reports_p, f'{output_tag}_final_report.html'),
                        f'{output_tag}_final_report',
                        runs_p=runs_p,
                        metadata_p=metadata_p,
                        geo=geo,
                        taxa=taxa)
            _set_tool_status(flags, 'final_report', 'completed')
        except Exception as exc:
            _set_tool_status(flags, 'final_report', 'failed', str(exc))
            flags['Tool Versions'] = dict(tool_versions)
            _write_flags_json(active_reports_p, output_tag, flags)
            raise
        flags['Tool Versions'] = dict(tool_versions)
        _write_flags_json(active_reports_p, output_tag, flags)

    if write_tool_versions_tsv and not batch and tool_versions:
        _write_tool_versions_tsv(active_reports_p, output_tag, tool_versions)

    print(f'Finished: {filename}')
    return {
        'flags': flags,
        'output_tag': output_tag,
        'mode': mode,
        'flumut_input_fp': flumut_input_fp if flags['Tool Status'].get('flumut', {}).get('status') == 'completed' else '',
        'batch_flumut_input_fp': batch_flumut_input_fp if flags['Tool Status'].get('flumut', {}).get('status') == 'completed' else '',
    }
def main(flagdict: dict = flagdict) -> None:
    """Main entry point for the pipeline."""
    args = parser()
    cwd = os.getcwd()
    project_root = Path(__file__).resolve().parent

    config_file = str(Path(args.config).expanduser().resolve())
    mode = args.mode
    update_flumut = args.update_flumut_db
    off_apps = args.turn_off
    force_apps = args.force
    single = args.single_sample.lower() == 'on'
    skip_cdhit = args.skip_cdhit
    requested_output_name = args.output_name
    write_tool_versions_tsv = args.tool_versions
    custom_output_root = resolve_output_root(cwd, args.outdir)
    replace_outputs = args.replace

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
    geo_json=config['Filenames']['geo']
    taxa_json=config['Filenames']['taxa']

    rm_previous = config['Functions']['remove_previous'] if args.remove_previous.lower() == 'on' else False

    samples_p = resolve_config_path(samples, project_root)
    runs_root_p = resolve_config_path(runs, project_root)
    references_p = resolve_config_path(references, project_root)
    reports_p = resolve_config_path(reports, project_root)
    logs_p = resolve_config_path(logs, project_root)
    blasts_p = resolve_config_path(blasts, project_root)
    clusters_p = resolve_config_path(clusters, project_root)
    metadata_p = resolve_config_path(metadata, project_root)

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
    print('Skip CD-HIT:', skip_cdhit)
    print('Batch directory:', args.batch_dir if args.batch_dir else '.')
    print('Output directory:', custom_output_root if custom_output_root else 'default config paths')
    print('Replace outputs:', replace_outputs)

    # Check required paths
    required_paths = [
        (samples_p, "Samples"),
        (runs_root_p, "Runs"),
        (references_p, "References"),
        (blasts_p, "Blast database"),
        (clusters_p, "Cluster database"),
        (metadata_p, "Metadata"),
    ]

    if args.batch or not custom_output_root:
        required_paths.append((reports_p, "Reports"))

    for pth, label in required_paths:
        if not os.path.exists(pth):
            print(f"{label} path does not exist: {pth}")
            sys.exit(1)

    if not os.path.exists(logs_p):
        os.makedirs(logs_p)

    # Optional cleanup
    if rm_previous:
        # Clean runs directory
        for file in glob.glob(os.path.join(runs_root_p, "*")):
            if os.path.isfile(file):
                os.remove(file)

        for run_dir in glob.glob(os.path.join(runs_root_p, "session_*")):
            if os.path.isdir(run_dir):
                shutil.rmtree(run_dir, ignore_errors=True)

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

    active_runs_p = create_run_workspace(runs_root_p)
    run_succeeded = False
    try:
        # Batch mode
        if args.batch:
            try:
                fasta_files = find_fasta_files(samples_p, batch_dir=args.batch_dir, recursive=True)
            except Exception as e:
                print(f"Could not scan batch directory: {e}")
                sys.exit(1)

            batch_label = args.batch_dir.strip().replace(os.sep, "_").replace("/", "_") if args.batch_dir else "batch_root"
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            try:
                batch_name = resolve_batch_name(args.batch_dir, timestamp, args.batch_name)
            except ValueError as exc:
                print(str(exc))
                sys.exit(1)

            exit_code = run_batch_pipeline(
                fasta_files=fasta_files,
                batch_dir_label=args.batch_dir,
                search_root_label=f"samples/{args.batch_dir}" if args.batch_dir else 'samples/.',
                batch_name=batch_name,
                config=config,
                batch_samples_p=samples_p,
                runs_p=active_runs_p,
                references_p=references_p,
                reports_p=reports_p,
                logs_p=logs_p,
                blasts_p=blasts_p,
                clusters_p=clusters_p,
                metadata_p=metadata_p,
                custom_output_root=custom_output_root,
                replace_outputs=replace_outputs,
                mode=mode,
                single=single,
                update_flumut=update_flumut,
                off_apps=off_apps,
                force_apps=force_apps,
                rm_previous=rm_previous,
                max_len=max_len,
                min_len=min_len,
                geo=geo_json,
                taxa=taxa_json,
                skip_cdhit=skip_cdhit,
                write_tool_versions_tsv=write_tool_versions_tsv,
            )
            if exit_code != 0:
                sys.exit(exit_code)
            run_succeeded = True
            return

        # Single-file mode
        try:
            resolved_filename = resolve_single_sample_file(samples_p, args.filename)
        except Exception as e:
            print(f"Could not resolve sample file: {e}")
            sys.exit(1)

        effective_output_name = requested_output_name
        if args.batch_fasta:
            temp_batch_dir = tempfile.mkdtemp(prefix='batch_fasta_', dir=active_runs_p)
            try:
                demux_report = demultiplex_batch_fasta(
                    os.path.join(samples_p, resolved_filename),
                    temp_batch_dir,
                    args.bf_regex,
                )

                if demux_report['demultipliable']:
                    temp_fasta_files = find_fasta_files(temp_batch_dir, recursive=True)
                    batch_label = make_sample_stem(resolved_filename)
                    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                    try:
                        batch_name = resolve_batch_name(batch_label, timestamp, args.batch_name)
                    except ValueError as exc:
                        print(str(exc))
                        sys.exit(1)

                    exit_code = run_batch_pipeline(
                        fasta_files=temp_fasta_files,
                        batch_dir_label=batch_label,
                        search_root_label='temporary batch-fasta workspace',
                        batch_name=batch_name,
                        config=config,
                        batch_samples_p=temp_batch_dir,
                        runs_p=active_runs_p,
                        references_p=references_p,
                        reports_p=reports_p,
                        logs_p=logs_p,
                        blasts_p=blasts_p,
                        clusters_p=clusters_p,
                        metadata_p=metadata_p,
                        custom_output_root=custom_output_root,
                        replace_outputs=replace_outputs,
                        mode=mode,
                        single=True,
                        update_flumut=update_flumut,
                        off_apps=off_apps,
                        force_apps=force_apps,
                        rm_previous=rm_previous,
                        max_len=max_len,
                        min_len=min_len,
                        geo=geo_json,
                        taxa=taxa_json,
                        skip_cdhit=skip_cdhit,
                        write_tool_versions_tsv=write_tool_versions_tsv,
                    )
                    if exit_code != 0:
                        sys.exit(exit_code)
                    run_succeeded = True
                    return

                if not args.force_bf:
                    unmatched = demux_report['unmatched_headers']
                    preview = ', '.join(unmatched[:3])
                    extra = '' if len(unmatched) <= 3 else f' (+{len(unmatched) - 3} more)'
                    raise ValueError(
                        f"Input FASTA is not demultipliable with --bf-regex; {len(unmatched)} header(s) failed to match: {preview}{extra}"
                    )

                if args.batch_name:
                    print('Ignoring --batch_name because demultiplexing failed and --force-bf switched to multi-sample mode.')
                effective_output_name = None
                single = False
            except Exception as e:
                shutil.rmtree(temp_batch_dir, ignore_errors=True)
                print(f"Batch FASTA preprocessing failed: {e}")
                sys.exit(1)
            finally:
                shutil.rmtree(temp_batch_dir, ignore_errors=True)

        sample_output_paths = resolve_single_sample_output_paths(
            resolved_filename,
            reports_p,
            references_p,
            custom_output_root,
            output_name=effective_output_name,
        )

        single_prepare_targets = [sample_output_paths['sample_root']] if custom_output_root else [
            sample_output_paths['sample_reports_p'],
            sample_output_paths['sample_reference_outputs_p'],
        ]

        try:
            prepare_output_targets(single_prepare_targets, replace_outputs)
        except FileExistsError as exc:
            print(str(exc))
            sys.exit(1)

        run_pipeline_for_file(
            filename=resolved_filename,
            flags_template=flagdict,
            config=config,
            samples_p=samples_p,
            runs_p=active_runs_p,
            references_p=references_p,
            reports_p=sample_output_paths['sample_reports_p'],
            logs_p=logs_p,
            blasts_p=blasts_p,
            clusters_p=clusters_p,
            metadata_p=metadata_p,
            reference_output_p=sample_output_paths['sample_reference_outputs_p'],
            mode=mode,
            single=single,
            update_flumut=update_flumut,
            off_apps=off_apps,
            force_apps=force_apps,
            rm_previous=rm_previous,
            max_len=max_len,
            min_len=min_len,
            geo=geo_json,
            taxa=taxa_json,
            skip_cdhit=skip_cdhit,
            batch=False,
            batch_dir='',
            timestamp='',
            output_name=effective_output_name,
            capture_tool_versions=True,
            write_tool_versions_tsv=write_tool_versions_tsv,
        )
        run_succeeded = True
    finally:
        if run_succeeded:
            shutil.rmtree(active_runs_p, ignore_errors=True)
if __name__=='__main__':
    main()
