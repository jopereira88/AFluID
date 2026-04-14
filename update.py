# AFluID update script
# Importing dependencies
import pandas as pd
import matplotlib.pyplot as plt
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
from main import fasta_preprocess, cd_hit_est_2d, cluster_assign, cluster_compile, cluster_miner, bclust, reblast, report_compiler


def update_ini_config(config_path, updates):
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
def normalise_segments(dataframe):
    '''
    Converts Segment values to int
    
    :param dataframe: A NCBI virus metadata dataframe (pandas DataFrame)

    Returns: dataframe
    '''
    dataframe['Segment']=dataframe['Segment'].astype(str)
    segment_to_number_iav={'PB2':'1', 'PB1':'2', 'PA':'3', 'HA':'4', 'NP':'5', 'NA':'6', 'MP':'7', 'NS':'8'}
    dataframe['Segment']=dataframe['Segment'].apply(lambda x: segment_to_number_iav[x] if x in segment_to_number_iav else x)
    dataframe['Segment']=dataframe['Segment'].astype(int)
    
    return dataframe

def remove_v_genotype(dataframe):
    '''
        Removes 'v' form genotype classifications
    
    :param dataframe: A NCBI virus metadata dataframe (pandas DataFrame)
    
    Returns: dataframe
    '''
    dataframe['Genotype']=dataframe['Genotype'].str.replace('v', '', regex=False)
    return dataframe

def normalise_len_to_int(dataframe):
    '''
    Converts Length values to int
    
    :param dataframe: A NCBI virus metadata dataframe (pandas DataFrame)

    Returns: dataframe
    '''
    dataframe['Length']=dataframe['Length'].astype(int)
    return dataframe


def normalize_to_list(value):
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


def extract_h_n_candidates(value_list):
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
    segment,
    subtype,
    segment_ha=4,
    segment_na=6,
    default_mode="keep"
):
    values = normalize_to_list(subtype)
    h_vals, n_vals = extract_h_n_candidates(values)

    if str(segment) == str(segment_ha):
        return h_vals[0] if h_vals else (subtype if default_mode == "keep" else None)

    if str(segment) == str(segment_na):
        return n_vals[0] if n_vals else (subtype if default_mode == "keep" else None)

    return subtype if default_mode == "keep" else None

def safe_str_to_dict(x):
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
def drop_na_metadata(dataframe):
    '''
    Drops na values from vital columns - Segment, Genotype, Host, Country
    
    :param dataframe: A NCBI virus metadata dataframe (pandas DataFrame)
    
    Returns: Dataframe without NA values for genotype, country, segment and host
    '''
    cols=['Segment','Genotype','Host','Country']
    dataframe= dataframe.dropna(subset=cols)
    dataframe= dataframe[(dataframe[cols] != "").all(axis=1)]
    return dataframe

def clean_assemblies(dataframe,number_of_segs):
    '''
    Removes every assembly that does not have the set number 
    of segments to be complete (default -IAV:8)
    
    :param dataframe: A NCBI virus metadata dataframe (pandas DataFrame)
    :param number_of_segs: Number of segments (integer)
    
    Returns: Dataframe
    '''
    assembly_counts = dataframe["Assembly"].value_counts()
    assemblies_to_remove = assembly_counts[assembly_counts < number_of_segs].index

    assembly_removal_mask = dataframe["Assembly"].isin(assemblies_to_remove)

    cleaned_df = dataframe.loc[~assembly_removal_mask].copy()
    return cleaned_df

def get_offending_accessions_genotype(df, assembly_col="Assembly", accession_col="Accession", genotype_col="Genotype"):
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

    def is_offending_genotype(g):
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

def get_size_outliers(dataframe, len_dict=seg_thress_A):
    '''
    Returns a list with size outlier sequences from a dataframe 
    
    :param dataframe: A NCBI virus metadata dataframe (pandas DataFrame)
    :param len_dict: A Dictionary with length thressholds (Dict)

    Returns: size_outliers (List)
    '''
    
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
    sequences_dict,
    ambiguous_nucleotides,
    threshold=10.0,
    header_delimiter="|"
):
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



def compare_non_determining(metadata_dataframe, report_dataframe):
    '''
    Compares metadata information to information obtained from an identification run 
    on AFluID.
    Tests if segments are consistent and if the genotype found on metadata is the 
    same as one of the genotypes of the sample background.
    Uses the accession number, parsed from the SAMPLE_NAME column of the report_dataframe 
    to join both df's for comparison.
    Assumes standard formats and cases (Upper for AFluID, casefold for NCBI / 
    >accession |metadata as sample name )
    This function is for non-determining segments only.

    :param metadata_dataframe: NCBI Virus metadata dataframe
    :param report_dataframe: AFluID ID report dataframe.

    Returns: List of inconsistent accessions.

    '''
    pos = report_dataframe.columns.get_loc('SAMPLE_NAME')+1
    report_dataframe.insert(pos, 'Accession', report_dataframe['SAMPLE_NAME'].apply(lambda x: x.split('|')[0].replace(' ', '').replace('>', '')))
    report_dataframe['SEGMENT']=report_dataframe['SEGMENT'].apply(safe_str_to_dict)
    report_dataframe['GENOTYPE']=report_dataframe['GENOTYPE'].apply(safe_str_to_dict)
    report_dataframe['SEGMENT']=report_dataframe['SEGMENT'].apply(lambda x: int(list(x.keys())[0]) if isinstance(x, dict) and len(x) == 1 else int(x))
    report_dataframe['GENOTYPE']=report_dataframe['GENOTYPE'].apply(lambda x: list(x.keys()) if isinstance(x, dict) else x)
    non_determining=report_dataframe[(report_dataframe['SEGMENT']!=4) & (report_dataframe['SEGMENT']!=6)]
    non_determining=non_determining.merge(metadata_dataframe, on='Accession', how='inner')

    offending_keys = []
    for key, v1, v2 in zip(non_determining["Accession"], non_determining["SEGMENT"], non_determining["Segment"]):
        if int(v1) != int(v2):
            offending_keys.append(key)
    for key, v1, v2 in zip(non_determining["Accession"], non_determining["GENOTYPE"], non_determining["Genotype"]):
        if v2 not in v2:
            offending_keys.append(key)
    
    return offending_keys

def compare_determining(metadata_dataframe, report_dataframe):
    '''
    Compares metadata information to information obtained from an identification run 
    on AFluID.
    Tests if segments are consistent and if the genotype found on metadata is the 
    same as one of the genotypes of the sample background.
    Uses the accession number, parsed from the SAMPLE_NAME column of the report_dataframe 
    to join both df's for comparison.
    Assumes standard formats and cases (Upper for AFluID, casefold for NCBI / 
    >accession |metadata as sample name )
    This function is for determining segments only.

    :param metadata_dataframe: NCBI Virus metadata dataframe
    :param report_dataframe: AFluID ID report dataframe.

    Returns: List of inconsistent accessions.

    '''
    pos = report_dataframe.columns.get_loc('SAMPLE_NAME')+1
    report_dataframe.insert(pos, 'Accession', report_dataframe['SAMPLE_NAME'].apply(lambda x: x.split('|')[0].replace(' ', '').replace('>', '')))
    report_dataframe['SEGMENT']=report_dataframe['SEGMENT'].apply(safe_str_to_dict)
    report_dataframe['GENOTYPE']=report_dataframe['GENOTYPE'].apply(safe_str_to_dict)
    report_dataframe['SEGMENT']=report_dataframe['SEGMENT'].apply(lambda x: int(list(x.keys())[0]) if isinstance(x, dict) and len(x) == 1 else int(x))
    report_dataframe['GENOTYPE']=report_dataframe['GENOTYPE'].apply(lambda x: list(x.keys()) if isinstance(x, dict) else x)
    determining=report_dataframe[(report_dataframe['SEGMENT']==4)|(report_dataframe['SEGMENT']==6)]
    determining=determining.merge(metadata_dataframe, on='Accession', how='inner')
    determining.loc[:,"GEN_EX"] = [expose_determining_component(seg,geno) for seg, geno in zip(determining["SEGMENT"],determining["GENOTYPE"])]
    determining.loc[:,"Gen_ex"] = [expose_determining_component(seg,geno) for seg, geno in zip(determining["Segment"],determining["Genotype"])]
    
    offending_keys=[]
    for key, v1, v2 in zip(determining["Accession"], determining["SEGMENT"], determining["Segment"]):
        if int(v1) != int(v2):
            offending_keys.append(key)
    for key, v1, v2 in zip(determining["Accession"], determining["Gen_ex"], determining["GEN_EX"]):
        if str(v1) != str(v2):
            offending_keys.append(key)
    
    return offending_keys

def remove_duplicates_list(in_list):
    '''
    Removes duplicate values from a list
    :param in_list: a list of values
    Returns dedup (list)

    '''
    dedup= set(in_list)
    return list(dedup)

def args():
    parser=argparse.ArgumentParser(prog='AFluID update routine',description='AFluID: Automated Influenza Identification Pipeline')
    parser.add_argument('-c','--config',type=str,help='Configuration file path',default='config.ini',required=False)
    parser.add_argument('-ff','--fastafile',type=str,help='Update fasta file name - path in config',required=True)
    parser.add_argument('-mf','--metadatafile',type=str,help='Update metadata file name - path in config',required=True)
    parser.add_argument('-rm','--remove_previous',type=str,help='Remove previous updates on update',choices=('on','off'),default='off',required=False)
    args=parser.parse_args()
    return args
    


def main():
    ##### LOADING VARIABLES
    arg=args()
    config_file=arg.config
    config=configparser.ConfigParser()
    config.read(config_file)
    flags=deepcopy(flagdict)
    main_metadata=config["Filenames"]["metadata"]
    fasta_filename=arg.fastafile
    old_fasta=f'{config["Filenames"]["l_blast"]}.fasta'
    new_metadata_filename=arg.metadatafile
    cwd=os.getcwd()
    update, runs, references, reports, logs = config['Paths']['update'], config['Paths']['runs'],\
          config['Paths']['references'], config['Paths']['reports'], config['Paths']['logs']
    blasts, clusters, metadata =config['Paths']['blast_database'], config['Paths']['cluster_database'],\
        config['Paths']['metadata']
    rm_previous=config['Functions']['remove_previous'] if\
        arg.remove_previous.lower()=='on' else False
    update_p, runs_p=os.path.abspath(os.path.join(cwd,update)), os.path.abspath(os.path.join(cwd,runs))
    references_p, reports_p=os.path.abspath(os.path.join(cwd,references)), os.path.abspath(os.path.join(cwd,reports))
    logs_p,blasts_p =os.path.abspath(os.path.join(cwd,logs)),os.path.abspath(os.path.join(cwd,blasts)) 
    clusters_p, metadata_p=os.path.abspath(os.path.join(cwd,clusters)), os.path.abspath(os.path.join(cwd,metadata))
    threads, num_seqs =int(config['blast']['num_threads']), int(config['blast']['max_target_seqs'])
    max,min=int(config["Sequence_Size"]["max"]), int(config["Sequence_Size"]["min"])
    date=datetime.now().strftime("%Y%m%d")
    ##### CHECKING PATHS
    if os.path.exists(update_p):
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
        print("Update path does not exist")
        sys.exit()
    if not os.path.exists(logs_p):
        os.makedirs(logs_p)
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

    #### STEP 4 - AFLUID RUN
    print("Starting AFluID")
    fasta_filename=fasta_filename+'.fasta'
    dict_cluster=json_load(os.path.join(clusters_p,config["Filenames"]["cluster_pkl"]))
    mappings=fasta_preprocess(fasta_filename,update_p,runs_p,min,max,flags,verbose=True)
    file=fasta_filename.replace('.fasta','')
    cd_hit_est_2d(os.path.join(runs_p,f'format_{fasta_filename}'),os.path.join(clusters_p,config["Filenames"]["cluster"]),\
                  os.path.join(runs_p,f'format_{file}'),float(config["CD-HIT"]["identity"]),logs_p,verbose=True)

    assignments=cluster_assign(f'format_{file}.clstr',os.path.join(runs_p,f'format_{fasta_filename}'),runs_p)
    cluster_compile(assignments,dict_cluster,os.path.join(runs_p,f'format_{file}.assign'),reports_p)
    cluster_report=cluster_miner(reports_p,file,runs_p,runs_p,flags)
    threads=config['blast']['num_threads']
    num_seqs=config['blast']['max_target_seqs']

    if flags['Master']['C_BLAST']:
        blast_cluster=bclust(clusters_p,config["Filenames"]["cluster_metadata"],runs_p,blasts_p,config["Filenames"]["cluster"],runs_p,\
                             fasta_filename,flags,num_threads=int(threads),max_tar_seq=int(num_seqs),fasta=f"{fasta_filename.replace('.fasta','')}_to_blast.fasta")
        if flags['Master']['L_BLAST']:
            blast_report=reblast(metadata_p,main_metadata,runs_p,blasts_p,config["Filenames"]["l_blast"],runs_p,\
                                 fasta_filename,threads=int(threads),max_tar_seq=int(num_seqs),fasta=f"{fasta_filename.replace('.fasta','')}_to_reblast.fasta")
            report_compiler(cluster_report,runs_p,fasta_filename,mappings,reports_p,flags,bclust_dict=blast_cluster,blast_dict=blast_report)
        else:
            report_compiler(cluster_report,runs_p,fasta_filename,mappings,reports_p,flags,bclust_dict=blast_cluster)
    else:
        report_compiler(cluster_report,runs_p,fasta_filename,mappings,reports_p,flags)
    
    #### STEP 5 - SEGMENT AND GENOTYPE COMPARISON
    print("Cross-referencing with metadata")
    id_report=pd.read_table(os.path.join(reports_p,f"{file}_ID_Report.txt"),index_col=False)
    accessions_to_remove=compare_non_determining(new_metadata,id_report)+compare_determining(new_metadata,id_report)
    accessions_to_remove=remove_duplicates_list(accessions_to_remove)
    new_metadata=new_metadata[~new_metadata['Accession'].isin(accessions_to_remove)]
    new_metadata=clean_assemblies(new_metadata,8)
    new_seqs={access_header[i]:update_sequences[access_header[i]] for i in new_metadata['Accession'].to_list() if i in access_header}

    #### STEP 6 - CONCATENATING FASTAS AND METADATA FILES
    ##### STEP 6.1 - BACKUP
    print("Creating backup and updating databases")
    backup_path='update_backup.tar.gz'
    source_files=[old_fasta,os.path.join(metadata_p,main_metadata),config_file]
    with tarfile.open(backup_path, "w:gz") as tar:
        for file_path in source_files:
            tar.add(file_path, arcname=os.path.basename(file_path))
    ##### STEP 6.2 - CONCATENATION
    with open(old_fasta, 'a') as f:
        for header, seq in new_seqs.items():
            f.write(f"{header}\n{seq}\n")
    os.rename(old_fasta,f'sequencesDnaInf_{date}'+'.fasta')
    rename_map={col:col.upper() for col in new_metadata.columns}
    new_metadata=new_metadata.rename(columns=rename_map)
    for col in main_metadata.columns:
        if col not in new_metadata.columns:
            new_metadata[col] = pd.NA
    new_metadata=new_metadata[main_metadata.columns]
    metadata_concat=pd.concat([main_metadata,new_metadata], ignore_index=True)
    metadata_concat.to_csv(f"flu_metadata_{date}_.csv",index=False,sep=';')

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
    os.system(f"python install.py -c {config_file}")

if __name__ == "__main__":
    main()

