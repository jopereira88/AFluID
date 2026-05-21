# AFluID update script
# Importing dependencies
import pandas as pd
from Bio import SeqIO, Entrez
from structures import seg_thress_A
import re
import ast
from flu_utils import seq_get, dict_to_fasta, json_load
from collections import Counter
import configparser
import argparse
from structures import flagdict,ambiguous_nucleotides
from copy import deepcopy, copy
import os
import sys
import glob
import tarfile
from datetime import datetime
from typing import Any, Iterable
from main import fasta_preprocess, cd_hit_est_2d, cluster_assign, cluster_compile, cluster_miner, bclust, reblast, report_compiler


def update_ini_config(config_path: str, updates: dict) -> None:
    """
    Update values in an INI configuration file.

    Parameters
    ----------
    config_path : str
        Path to the INI file.
    updates : dict
        Nested dictionary in the form:
        {
            "section": {
                "key": "value"
            }
        }

    Notes
    -----
    Missing sections are created automatically.
    All values are written as strings.
    """
    config = configparser.ConfigParser()
    config.read(config_path, encoding="utf-8")

    for section, section_updates in updates.items():
        if not config.has_section(section):
            config.add_section(section)

        for key, value in section_updates.items():
            config[section][key] = str(value)

    with open(config_path, "w", encoding="utf-8") as f:
        config.write(f)


# Value normalization and transformation funcs
def normalise_segments(dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Convert segment labels to integer segment numbers.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        NCBI Virus metadata table containing a ``Segment`` column.

    Returns
    -------
    pandas.DataFrame
        Input dataframe with normalized integer values in ``Segment``.
    """
    dataframe['Segment']=dataframe['Segment'].astype(str)
    segment_to_number_iav={'PB2':'1', 'PB1':'2', 'PA':'3', 'HA':'4', 'NP':'5', 'NA':'6', 'MP':'7', 'NS':'8'}
    dataframe['Segment']=dataframe['Segment'].apply(lambda x: segment_to_number_iav[x] if x in segment_to_number_iav else x)
    dataframe['Segment']=dataframe['Segment'].astype(int)
    
    return dataframe

def remove_v_genotype(dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Remove trailing ``v`` markers from genotype values.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        NCBI Virus metadata table containing a ``Genotype`` column.

    Returns
    -------
    pandas.DataFrame
        Input dataframe with normalized genotype strings.
    """
    dataframe['Genotype']=dataframe['Genotype'].str.replace('v', '', regex=False)
    return dataframe

def normalise_len_to_int(dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Convert sequence lengths to integers.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        NCBI Virus metadata table containing a ``Length`` column.

    Returns
    -------
    pandas.DataFrame
        Input dataframe with integer values in ``Length``.
    """
    dataframe['Length']=dataframe['Length'].astype(int)
    return dataframe


def normalize_to_list(value: Any) -> list[Any]:
    """
    Normalize a scalar or serialized list into a Python list.

    Parameters
    ----------
    value : Any
        Value to normalize.

    Returns
    -------
    list
        Empty list for null values, the original list for list inputs, or a
        single-item list for scalar values.
    """
    if value is None:
        return []

    if isinstance(value, float) and pd.isna(value):
        return []

    if isinstance(value, list):
        return value

    if isinstance(value, str):
        value = value.strip()

        try:
            parsed = ast.literal_eval(value)
            if isinstance(parsed, list):
                return parsed
        except (ValueError, SyntaxError):
            pass

        return [value]

    return [value]


def extract_h_n_candidates(value_list: Iterable[Any]) -> tuple[list[str], list[str]]:
    """
    Extract H and N subtype tokens from a list of labels.

    Parameters
    ----------
    value_list : list
        Candidate genotype strings such as ``H5N1``, ``H5``, or ``N1``.

    Returns
    -------
    tuple[list, list]
        Lists containing the H and N subtype tokens found in the input.
    """
    h_vals = []
    n_vals = []

    for item in value_list:
        if not isinstance(item, str):
            continue

        item = item.strip().upper()

        # normalizar Ñ -> N
        item = item.replace("Ñ", "N")

        # remover 'V' final em casos como H1N2V
        item = re.sub(r'V$', '', item)

        # caso completo: H1N2 ou H1N2V
        m = re.fullmatch(r'(H\d+)(N\d+)', item)
        if m:
            h_vals.append(m.group(1))
            n_vals.append(m.group(2))
            continue

        # caso isolado: H1
        if re.fullmatch(r'H\d+', item):
            h_vals.append(item)
            continue

        # caso isolado: N2
        if re.fullmatch(r'N\d+', item):
            n_vals.append(item)
            continue

    return h_vals, n_vals


def expose_determining_component(
    segment: Any,
    subtype: Any,
    segment_ha: int = 4,
    segment_na: int = 6,
    default_mode: str = "keep"
) -> Any:
    """
    Extract the determining H or N component for a segment.

    Parameters
    ----------
    segment : Any
        Segment number being evaluated.
    subtype : Any
        Genotype value or collection of genotype values.
    segment_ha : int, default 4
        Segment number representing HA.
    segment_na : int, default 6
        Segment number representing NA.
    default_mode : {"keep", "drop"}, default "keep"
        Fallback behavior when a determining component cannot be extracted.

    Returns
    -------
    str or None
        The extracted H or N subtype token, the original value, or None
        depending on the segment and fallback mode.
    """
    values = normalize_to_list(subtype)
    h_vals, n_vals = extract_h_n_candidates(values)

    if str(segment) == str(segment_ha):
        return h_vals[0] if h_vals else (subtype if default_mode == "keep" else None)

    if str(segment) == str(segment_na):
        return n_vals[0] if n_vals else (subtype if default_mode == "keep" else None)

    return subtype if default_mode == "keep" else None

def safe_str_to_dict(x: Any) -> Any:
    """
    Parse a string into a dictionary only when it represents one.

    Parameters
    ----------
    x : Any
        Value to inspect.

    Returns
    -------
    Any
        Parsed dictionary when ``x`` is a string representation of a dict,
        otherwise the original value.
    """
    # does not eval dicts
    if isinstance(x, dict):
        return x
    
    # will only convert strings
    if not isinstance(x, str):
        return x
    
    try:
        value = ast.literal_eval(x)
        return value if isinstance(value, dict) else x
    except (ValueError, SyntaxError):
        return x
    
# Filtering funcs
def drop_na_metadata(dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Remove rows missing required metadata fields.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        NCBI Virus metadata table.

    Returns
    -------
    pandas.DataFrame
        Filtered dataframe with non-null and non-empty values in the required
        columns.
    """
    cols=['Segment','Genotype','Host','Country']
    dataframe= dataframe.dropna(subset=cols)
    dataframe= dataframe[(dataframe[cols] != "").all(axis=1)]
    return dataframe

def clean_assemblies(dataframe: pd.DataFrame,number_of_segs: int) -> pd.DataFrame:
    """
    Keep only assemblies with the expected number of segments.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        NCBI Virus metadata table containing an ``Assembly`` column.
    number_of_segs : int
        Minimum number of segments required for an assembly to be retained.

    Returns
    -------
    pandas.DataFrame
        Filtered dataframe containing only complete assemblies.
    """
    assembly_counts = dataframe["Assembly"].value_counts()
    assemblies_to_remove = assembly_counts[assembly_counts < number_of_segs].index

    assembly_removal_mask = dataframe["Assembly"].isin(assemblies_to_remove)

    cleaned_df = dataframe.loc[~assembly_removal_mask].copy()
    return cleaned_df

def get_offending_accessions_genotype(df: pd.DataFrame, assembly_col: str = "Assembly", accession_col: str = "Accession", genotype_col: str = "Genotype") -> list[str]:
    """
    Returns a list of accessions corresponding to assemblies with incomplete, mixed, empty, or non-standard genotypes.

    Regra:
    - An assembly offends if:
      - incomplete (ex.: H5, N1)
      - mixed
      - empty / NaN
      - any value not specifically HxNy
    - whenever an assembly has ANY offending genotype,
      the whole assembly is considered offending, and all its accessions are returned
    - The original df will not be modified, and the returned list will contain unique accessions without NaN values.
    """

    def is_offending_genotype(g: Any) -> bool:
        if pd.isna(g):
            return True

        s = str(g).strip().casefold()

        if s == "":
            return True

        if "mixed" in s:
            return True

        # Only considering the full HxNy format
        return re.fullmatch(r"h\d+n\d+", s) is None

    offending_assemblies = (
        df.groupby(assembly_col)[genotype_col]
          .apply(lambda genotypes: genotypes.apply(is_offending_genotype).any())
          .pipe(lambda s: s[s].index)
    )

    offending_accessions = (
        df.loc[df[assembly_col].isin(offending_assemblies), accession_col]
          .dropna()
          .drop_duplicates()
          .tolist()
    )

    return offending_accessions

def get_size_outliers(dataframe: pd.DataFrame, len_dict: dict = seg_thress_A) -> list[str]:
    """
    Return accessions whose segment lengths fall outside expected ranges.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        NCBI Virus metadata table containing ``Segment``, ``Length``, and
        ``Accession`` columns.
    len_dict : dict, default ``seg_thress_A``
        Length threshold dictionary keyed by segment-specific min/max names.

    Returns
    -------
    list
        Accessions whose segment length is outside the configured thresholds.
    """
    
    size_outliers_PA=dataframe[(dataframe['Segment']==3) & ((dataframe['Length']<len_dict['PA_min']) | (dataframe['Length']>len_dict['PA_max']))]
    size_outliers_PA=size_outliers_PA['Accession'].tolist()
    size_outliers_PB2=dataframe[(dataframe['Segment']==1) & ((dataframe['Length']<len_dict['PB2_min']) | (dataframe['Length']>len_dict['PB2_max']))]
    size_outliers_PB2=size_outliers_PB2['Accession'].tolist()
    size_outliers_PB1=dataframe[(dataframe['Segment']==2) & ((dataframe['Length']<len_dict['PB1_min']) | (dataframe['Length']>len_dict['PB1_max']))]
    size_outliers_PB1=size_outliers_PB1['Accession'].tolist()
    size_outliers_HA=dataframe[(dataframe['Segment']==4) & ((dataframe['Length']<len_dict['HA_min']) | (dataframe['Length']>len_dict['HA_max']))]
    size_outliers_HA=size_outliers_HA['Accession'].tolist()
    size_outliers_NP=dataframe[(dataframe['Segment']==5) & ((dataframe['Length']<len_dict['NP_min']) | (dataframe['Length']>len_dict['NP_max']))]
    size_outliers_NP=size_outliers_NP['Accession'].tolist()
    size_outliers_NA=dataframe[(dataframe['Segment']==6) & ((dataframe['Length']<len_dict['NA_min']) | (dataframe['Length']>len_dict['NA_max']))]
    size_outliers_NA=size_outliers_NA['Accession'].tolist()
    size_outliers_MP=dataframe[(dataframe['Segment']==7) & ((dataframe['Length']<len_dict['MP_min']) | (dataframe['Length']>len_dict['MP_max']))]
    size_outliers_MP=size_outliers_MP['Accession'].tolist()
    size_outliers_NS=dataframe[(dataframe['Segment']==8) & ((dataframe['Length']<len_dict['NS_min']) | (dataframe['Length']>len_dict['NS_max']))]
    size_outliers_NS=size_outliers_NS['Accession'].tolist()
    
    size_outliers=size_outliers_HA+size_outliers_MP+size_outliers_NA+size_outliers_NP+size_outliers_NS+size_outliers_PA+size_outliers_PB1+size_outliers_PB2
    return size_outliers

def filter_by_ambiguous_nucleotides(
    sequences_dict: dict[str, str],
    ambiguous_nucleotides: dict | set,
    threshold: float = 10.0,
    header_delimiter: str = "|"
) -> dict[str, float]:
    """
    Identify sequences whose percentage of ambiguous nucleotides exceeds a given threshold.

    Parameters
    ----------
    sequences_dict : dict
        Dictionary mapping sequence headers to nucleotide sequences.
        Example:
        {
            ">seq1|meta": "ATGCRYN...",
            ">seq2|meta": "ATGC..."
        }

    ambiguous_nucleotides : dict or set
        Collection of ambiguous nucleotide symbols.
        If a dictionary is provided, its keys are used.

    threshold : float, default 10.0
        Minimum percentage of ambiguous nucleotides required for a sequence
        to be included in the output.

    header_delimiter : str, default "|"
        Delimiter used to split the sequence header in order to extract the
        sequence identifier.

    Returns
    -------
    dict
        Dictionary mapping cleaned sequence identifiers to the percentage of
        ambiguous nucleotides for sequences exceeding the threshold.

    Notes
    -----
    The sequence identifier is extracted from the header by:
        1. splitting on `header_delimiter`,
        2. keeping the first field,
        3. removing ">" and spaces.

    Percentages are rounded to two decimal places before comparison.
    """
    ambiguous_set = (
        set(ambiguous_nucleotides.keys())
        if isinstance(ambiguous_nucleotides, dict)
        else set(ambiguous_nucleotides)
    )

    perc_ambiguous = {}

    for header, sequence in sequences_dict.items():
        if not sequence:
            continue

        counts = Counter(sequence)
        ambiguous_count = sum(counts[base] for base in ambiguous_set if base in counts)
        ambiguous_percent = round((ambiguous_count / len(sequence)) * 100, 2)

        if ambiguous_percent > threshold:
            seq_id = header.split(header_delimiter, 1)[0].replace(">", "").replace(" ", "")
            perc_ambiguous[seq_id] = ambiguous_percent

    return perc_ambiguous



def compare_non_determining(metadata_dataframe: pd.DataFrame, report_dataframe: pd.DataFrame) -> list[str]:
    """
    Compare report assignments against metadata for non-determining segments.

    Parameters
    ----------
    metadata_dataframe : pandas.DataFrame
        Reference metadata table from NCBI Virus.
    report_dataframe : pandas.DataFrame
        AFluID ID report dataframe.

    Returns
    -------
    list
        Accessions with inconsistent segment or genotype assignments.
    """
    report_dataframe = report_dataframe.copy()

    report_dataframe['Accession'] = report_dataframe['SAMPLE_NAME'].apply(
        lambda x: str(x).split('|')[0].replace(' ', '').replace('>', '')
    )

    report_dataframe['SEGMENT'] = report_dataframe['SEGMENT'].apply(safe_str_to_dict)
    report_dataframe['GENOTYPE'] = report_dataframe['GENOTYPE'].apply(safe_str_to_dict)

    def parse_segment(x: Any) -> int:
        if isinstance(x, dict) and len(x) == 1:
            return int(list(x.keys())[0])
        return int(x)

    def parse_genotype(x: Any) -> list[str]:
        if isinstance(x, dict):
            return [str(k) for k in x.keys()]
        elif isinstance(x, list):
            return [str(v) for v in x]
        elif pd.isna(x):
            return []
        else:
            return [str(x)]

    report_dataframe['SEGMENT'] = report_dataframe['SEGMENT'].apply(parse_segment)
    report_dataframe['GENOTYPE'] = report_dataframe['GENOTYPE'].apply(parse_genotype)

    non_determining = report_dataframe[
        (report_dataframe['SEGMENT'] != 4) & (report_dataframe['SEGMENT'] != 6)
    ].copy()

    non_determining = non_determining.merge(metadata_dataframe, on='Accession', how='inner')

    offending_keys = set()

    for key, v1, v2 in zip(non_determining["Accession"], non_determining["SEGMENT"], non_determining["Segment"]):
        if int(v1) != int(v2):
            offending_keys.add(key)

    for key, v1, v2 in zip(non_determining["Accession"], non_determining["GENOTYPE"], non_determining["Genotype"]):
        if str(v2) not in [str(g) for g in v1]:
            offending_keys.add(key)

    return list(offending_keys)


def compare_determining(metadata_dataframe: pd.DataFrame, report_dataframe: pd.DataFrame) -> list[str]:
    """
    Compare report assignments against metadata for determining segments.

    Parameters
    ----------
    metadata_dataframe : pandas.DataFrame
        Reference metadata table from NCBI Virus.
    report_dataframe : pandas.DataFrame
        AFluID ID report dataframe.

    Returns
    -------
    list
        Accessions with inconsistent determining segment assignments.
    """
    report_dataframe = report_dataframe.copy()

    report_dataframe['Accession'] = report_dataframe['SAMPLE_NAME'].apply(
        lambda x: str(x).split('|')[0].replace(' ', '').replace('>', '')
    )

    report_dataframe['SEGMENT'] = report_dataframe['SEGMENT'].apply(safe_str_to_dict)
    report_dataframe['GENOTYPE'] = report_dataframe['GENOTYPE'].apply(safe_str_to_dict)

    def parse_segment(x: Any) -> int:
        if isinstance(x, dict) and len(x) == 1:
            return int(list(x.keys())[0])
        return int(x)

    def parse_genotype(x: Any) -> list[str]:
        if isinstance(x, dict):
            return [str(k) for k in x.keys()]
        elif isinstance(x, list):
            return [str(v) for v in x]
        elif pd.isna(x):
            return []
        else:
            return [str(x)]

    report_dataframe['SEGMENT'] = report_dataframe['SEGMENT'].apply(parse_segment)
    report_dataframe['GENOTYPE'] = report_dataframe['GENOTYPE'].apply(parse_genotype)

    determining = report_dataframe[
        (report_dataframe['SEGMENT'] == 4) | (report_dataframe['SEGMENT'] == 6)
    ].copy()

    determining = determining.merge(metadata_dataframe, on='Accession', how='inner')

    determining.loc[:, "GEN_EX"] = [
        expose_determining_component(seg, geno)
        for seg, geno in zip(determining["SEGMENT"], determining["GENOTYPE"])
    ]
    determining.loc[:, "Gen_ex"] = [
        expose_determining_component(seg, geno)
        for seg, geno in zip(determining["Segment"], determining["Genotype"])
    ]

    offending_keys = set()

    for key, v1, v2 in zip(determining["Accession"], determining["SEGMENT"], determining["Segment"]):
        if int(v1) != int(v2):
            offending_keys.add(key)

    for key, v1, v2 in zip(determining["Accession"], determining["Gen_ex"], determining["GEN_EX"]):
        if str(v1) != str(v2):
            offending_keys.add(key)

    return list(offending_keys)

def remove_duplicates_list(in_list: list[Any]) -> list[Any]:
    """
    Remove duplicate values from a list.

    Parameters
    ----------
    in_list : list
        Values to deduplicate.

    Returns
    -------
    list
        Deduplicated values.
    """
    dedup= set(in_list)
    return list(dedup)

def args() -> argparse.Namespace:
    """
    Parse command-line arguments for the update routine.

    Returns
    -------
    argparse.Namespace
        Parsed update-script command-line options.
    """
    parser=argparse.ArgumentParser(prog='AFluID update routine',description='AFluID: Automated Influenza Identification Pipeline')
    parser.add_argument('-c','--config',type=str,help='Configuration file path',default='config.ini',required=False)
    parser.add_argument('-ff','--fastafile',type=str,help='Update fasta file name - path in config',required=True)
    parser.add_argument('-mf','--metadatafile',type=str,help='Update metadata file name - path in config',required=True)
    parser.add_argument('-rm','--remove_previous',type=str,help='Remove previous updates on update',choices=('on','off'),default='off',required=False)
    args=parser.parse_args()
    return args
    


def main() -> None:
    """
    Run the database update workflow.

    Returns
    -------
    None

    Notes
    -----
    This routine cleans incoming metadata and sequences, validates them through
    an AFluID run, merges accepted records into the local databases, updates the
    configuration, and re-runs the installation routine.
    """
    ##### LOADING VARIABLES
    arg=args()
    config_file=arg.config
    config=configparser.ConfigParser()
    config.read(config_file)
    flags=deepcopy(flagdict)
    fasta_filename=arg.fastafile
    old_fasta=f'{config["Filenames"]["l_blast"]}.fasta'
    new_metadata_filename=arg.metadatafile
    cwd=os.getcwd()
    if 'update_reports' not in config['Paths']:
        config['Paths']['update_reports'] = 'update/reports'
    update, runs, references, update_reports, logs = config['Paths']['update'], config['Paths']['runs'],\
          config['Paths']['references'], config['Paths']['update_reports'], config['Paths']['logs']
    blasts, clusters, metadata =config['Paths']['blast_database'], config['Paths']['cluster_database'],\
        config['Paths']['metadata']
    rm_previous=config['Functions']['remove_previous'] if\
        arg.remove_previous.lower()=='on' else False
    update_p, runs_p=os.path.abspath(os.path.join(cwd,update)), os.path.abspath(os.path.join(cwd,runs))
    references_p = os.path.abspath(os.path.join(cwd,references))
    update_reports_p = os.path.abspath(os.path.join(cwd, update_reports))
    logs_p,blasts_p =os.path.abspath(os.path.join(cwd,logs)),os.path.abspath(os.path.join(cwd,blasts)) 
    clusters_p, metadata_p=os.path.abspath(os.path.join(cwd,clusters)), os.path.abspath(os.path.join(cwd,metadata))
    threads, num_seqs =int(config['blast']['num_threads']), int(config['blast']['max_target_seqs'])
    max,min=int(config["Sequence_Size"]["max"]), int(config["Sequence_Size"]["min"])
    date=datetime.now().strftime("%Y%m%d")
    main_metadata=pd.read_csv(os.path.join(metadata_p,config["Filenames"]["metadata"]), sep=';', index_col=False)
    ##### CHECKING PATHS
    if os.path.exists(update_p):
        if os.path.exists(runs_p):
            if os.path.exists(references_p):
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
                print("References path does not exist")
                sys.exit()
        else:
            print("Runs path does not exist")
            sys.exit()
    else:
        print("Update path does not exist")
        sys.exit()
    if not os.path.exists(logs_p):
        os.makedirs(logs_p)
    if not os.path.exists(update_reports_p):
        os.makedirs(update_reports_p)
    if rm_previous:
        # Remove all files in RUN_DIR except *.pkl
        for file in glob.glob(os.path.join(runs_p, "*")):
            os.remove(file)
    #### STEP 1 - CLEANING UPDATE METADATA
    print('Cleaning Metadata')
    new_metadata=pd.read_csv(os.path.join(update_p,new_metadata_filename),index_col=False,low_memory=False)
    new_metadata=drop_na_metadata(new_metadata)
    new_metadata=normalise_segments(new_metadata)
    new_metadata=remove_v_genotype(new_metadata)
    new_metadata=normalise_len_to_int(new_metadata)

    #### STEP 2 - FILTERING BY SIZE AND GENOTYPE
    print('Filtering accessisons')
    accessions_to_remove=get_offending_accessions_genotype(new_metadata) + get_size_outliers(new_metadata)
    accessions_to_remove=remove_duplicates_list(accessions_to_remove)
    new_metadata=new_metadata[~new_metadata['Accession'].isin(accessions_to_remove)]
    new_metadata=clean_assemblies(new_metadata,8)

    #### STEP 3 - FILTERING BY AMBIGUOUS NUCLEOTIDE CONTENT AND FAST OUTPUT
    print('Checking sequence quality')
    update_sequences=seq_get(os.path.join(update_p, fasta_filename))
    access_header={header.split('|')[0].replace('>','').replace(' ','') : header for header in update_sequences}
    ambiguous=filter_by_ambiguous_nucleotides(update_sequences,ambiguous_nucleotides)
    ambiguous=list(ambiguous.keys())
    new_metadata=new_metadata[~new_metadata['Accession'].isin(ambiguous)]
    new_metadata=clean_assemblies(new_metadata,8)
    new_seqs={access_header[i]:update_sequences[access_header[i]] for i in new_metadata['Accession'].to_list() if i in access_header}
    fasta_filename=fasta_filename.replace('.fasta','_new')
    dict_to_fasta(new_seqs,os.path.join(update_p,fasta_filename))
    file=fasta_filename.replace('.fasta','')
    #### STEP 4 - AFLUID RUN
    if not os.path.exists(os.path.join(update_reports_p,f"{file}_ID_Report.txt")):
        print("Starting AFluID")
        fasta_filename=fasta_filename+'.fasta'
        dict_cluster=json_load(os.path.join(clusters_p,config["Filenames"]["cluster_pkl"]))
        print("Remapping seqs")
        mappings=fasta_preprocess(fasta_filename,update_p,runs_p,file,min,max,flags,verbose=True)

        #snapshot 1 for the longest procecesses
        if not os.path.exists(os.path.join(runs_p,f'format_{file}')):
            print("Running cd-hit")
            cd_hit_est_2d(os.path.join(runs_p,f'format_{fasta_filename}'),os.path.join(clusters_p,config["Filenames"]["cluster"]),\
                  os.path.join(runs_p,f'format_{file}'),float(config["CD-HIT"]["identity"]),logs_p,verbose=True)
        print('Creating assignment object')
        assignments=cluster_assign(f'format_{file}.clstr',os.path.join(runs_p,f'format_{fasta_filename}'),runs_p)
        print("Creating cluster report")
        cluster_compile(assignments,dict_cluster,os.path.join(runs_p,f'format_{file}.assign'),update_reports_p)
        print("Creating Blast Fasta")
        cluster_report=cluster_miner(update_reports_p,file,runs_p,runs_p,flags)
        threads=config['blast']['num_threads']
        num_seqs=config['blast']['max_target_seqs']
        print("Starting Blast")
        if flags['Master']['C_BLAST']:
            blast_cluster=bclust(clusters_p,config["Filenames"]["cluster_metadata"],runs_p,blasts_p,config["Filenames"]["cluster"],runs_p,\
                             file,flags,num_threads=int(threads),max_tar_seq=int(num_seqs),fasta=f"{file}_to_blast.fasta")
            if flags['Master']['L_BLAST']:
                blast_report=reblast(metadata_p,config["Filenames"]["metadata"],runs_p,blasts_p,config["Filenames"]["l_blast"],runs_p,\
                                 file,threads=int(threads),max_tar_seq=int(num_seqs),fasta=f"{file}_to_reblast.fasta")
                report_compiler(cluster_report,runs_p,file,mappings,update_reports_p,flags,bclust_dict=blast_cluster,blast_dict=blast_report)
            else:
                report_compiler(cluster_report,runs_p,file,mappings,update_reports_p,flags,bclust_dict=blast_cluster)
        else:
            report_compiler(cluster_report,runs_p,file,mappings,update_reports_p,flags)
    
    #### STEP 5 - SEGMENT AND GENOTYPE COMPARISON
    print("Cross-referencing with metadata")
    id_report=pd.read_table(os.path.join(update_reports_p,f"{file}_ID_Report.txt"),index_col=False)
    accessions_to_remove=compare_non_determining(new_metadata,id_report)+compare_determining(new_metadata,id_report)
    accessions_to_remove=remove_duplicates_list(accessions_to_remove)
    new_metadata=new_metadata[~new_metadata['Accession'].isin(accessions_to_remove)]
    new_metadata=clean_assemblies(new_metadata,8)
    new_seqs={access_header[i]:update_sequences[access_header[i]] for i in new_metadata['Accession'].to_list() if i in access_header}

    #### STEP 6 - CONCATENATING FASTAS AND METADATA FILES
    ##### STEP 6.1 - BACKUP
    if os.listdir(metadata_p) and not os.path.exists(f'flu_metadata_{date}_.csv'):
        print("Creating backup and updating databases")
        backup_path='update_backup.tar.gz'
        source_files=[old_fasta,os.path.join(metadata_p,config["Filenames"]["metadata"]),config_file]
        with tarfile.open(backup_path, "w:gz") as tar:
            for file_path in source_files:
                tar.add(file_path, arcname=os.path.basename(file_path))
    ##### STEP 6.2 - CONCATENATION
        with open(old_fasta, 'a') as f:
            for header, seq in new_seqs.items():
                f.write(f"{header}\n{seq}\n")
        rename_map={col:col.upper() for col in new_metadata.columns}
        new_metadata=new_metadata.rename(columns=rename_map)
        for col in main_metadata.columns:
            if col not in new_metadata.columns:
                new_metadata[col] = pd.NA
        new_metadata=new_metadata[main_metadata.columns]
        metadata_concat=pd.concat([main_metadata,new_metadata], ignore_index=True)
        metadata_concat.to_csv(f"flu_metadata_{date}_.csv",index=False,sep=';')
        os.rename(old_fasta,f'sequencesDnaInf_{date}'+'.fasta')
    #### STEP 7 - INSTALLATION
    ##### STEP 7.1 - UPDATING FILENAMES
        print("Regenerating databases")
        config_updates = {
        "Filenames": {
            "l_blast": f"sequencesDnaInf_{date}",
            "metadata": f"flu_metadata_{date}_.csv"
        },
        "Versions": {
            "update_date": date
        }
        }
        for file in glob.glob(os.path.join(blasts_p, "*")):
            os.remove(file)
        for file in glob.glob(os.path.join(clusters_p, "*")):
            os.remove(file)
        for file in glob.glob(os.path.join(metadata_p, "*")):
            os.remove(file)
        update_ini_config(config_file, config_updates)
    
    ###### STEP 7.2 - RE-RUNNING INSTALL ROUTINE
    os.system(f"python3 install.py {config_file}")

if __name__ == "__main__":
    main()
