import os
import sys
import configparser
import shutil
import subprocess
import csv
import re
from datetime import date
from flu_utils import parse_clstr
from collections import Counter
import json
from pathlib import Path


def _normalize_header(name: str) -> str:
    return name.strip().upper()


def _resolve_config_path(raw_path: str, project_root: Path) -> str:
    """Resolve config paths relative to the project root unless already absolute."""
    path = Path(raw_path).expanduser()
    if path.is_absolute():
        return str(path.resolve())
    return str((project_root / path).resolve())


def _parse_collection_date(value: str) -> tuple[date, str] | None:
    """Return a sortable date plus the original precision-preserving string."""
    if value is None:
        return None

    raw = str(value).strip()
    if not raw:
        return None

    match = re.fullmatch(r'(\d{4})(?:-(\d{2})(?:-(\d{2}))?)?', raw)
    if not match:
        return None

    year = int(match.group(1))
    month = int(match.group(2)) if match.group(2) else 1
    day = int(match.group(3)) if match.group(3) else 1

    try:
        normalized = date(year, month, day)
    except ValueError:
        return None

    return normalized, raw


def _format_collection_date_range(values: list[str]) -> str:
    parsed_values = []
    for value in values:
        parsed = _parse_collection_date(value)
        if parsed is not None:
            parsed_values.append(parsed)

    if not parsed_values:
        return 'Unknown'

    parsed_values.sort(key=lambda item: item[0])
    min_display = parsed_values[0][1]
    max_display = parsed_values[-1][1]

    if min_display == max_display:
        return min_display

    return f'{min_display} - {max_display}'


def _load_cluster_metadata_rows(metadata_fp: str, sep: str = ';') -> dict[str, list[str]]:
    """Load only the metadata fields needed for cluster aggregation."""
    required_fields = ['GENOTYPE', 'SEGMENT', 'HOST', 'COUNTRY']
    rows = {}

    with open(metadata_fp, 'r', newline='') as metadata_file:
        reader = csv.DictReader(metadata_file, delimiter=sep)
        if reader.fieldnames is None:
            raise ValueError(f'Metadata file {metadata_fp} is missing a header row.')

        header_map = {_normalize_header(name): name for name in reader.fieldnames}
        missing = [field for field in required_fields if field not in header_map]
        if missing:
            raise ValueError(
                f'Metadata file {metadata_fp} is missing required columns: {", ".join(missing)}'
            )

        accession_key = reader.fieldnames[0]
        for row in reader:
            accession = str(row.get(accession_key, '')).strip()
            if not accession:
                continue
            rows[accession] = [
                str(row.get(header_map['GENOTYPE'], '')).strip(),
                str(row.get(header_map['SEGMENT'], '')).strip(),
                str(row.get(header_map['HOST'], '')).strip(),
                str(row.get(header_map['COUNTRY'], '')).strip(),
                str(row.get(header_map['COLLECTION_DATE'], '')).strip() if 'COLLECTION_DATE' in header_map else '',
            ]

    return rows

def cluster_charact(metadata: str,metadata_p: str,clstr: str,cluster_p: str,cluster_pkl_name: str,cluster_meta_name: str) -> None:
    """
    Build cluster-level metadata summaries from a CD-HIT cluster file.

    Parameters
    ----------
    metadata : str
        Metadata filename to read from ``metadata_p``.
    metadata_p : str
        Directory containing the metadata file.
    clstr : str
        CD-HIT ``.clstr`` filename to read from ``cluster_p``.
    cluster_p : str
        Directory containing cluster artifacts and output files.
    cluster_pkl_name : str
        Stem used for the JSON summary output.
    cluster_meta_name : str
        Filename for the tabular cluster metadata report.

    Returns
    -------
    None

    Notes
    -----
    This function writes both a JSON summary and a tab-delimited cluster
    metadata table.
    """
    meta = _load_cluster_metadata_rows(os.path.join(metadata_p, metadata), sep=';')
    clusters=parse_clstr(os.path.join(cluster_p,clstr))
    #getting metadata annotations per cluster
    cluster_table={}
    for key in clusters:
        cluster_table[key]=[Counter(),Counter(),Counter(),Counter(),[]] #genotype,segment,host,country,collection_date
        for value in clusters[key]:
            cluster_table[key][0][meta[value][0]]+=1
            cluster_table[key][1][meta[value][1]]+=1
            cluster_table[key][2][meta[value][2]]+=1
            cluster_table[key][3][meta[value][3]]+=1
            cluster_table[key][4].append(meta[value][4])
    for elem in cluster_table:
        cluster_table[elem][0]=dict(cluster_table[elem][0])
        cluster_table[elem][1]=dict(cluster_table[elem][1])
        cluster_table[elem][2]=dict(cluster_table[elem][2])
        cluster_table[elem][3]=dict(cluster_table[elem][3])
        cluster_table[elem][4]=_format_collection_date_range(cluster_table[elem][4])
    #getting cluster representatives
    clusters=parse_clstr(os.path.join(cluster_p,clstr),access_only=False)
    cluster_reps={}
    for cl in clusters:
        rep=[item for item in clusters[cl] if '*' in item]
        cluster_reps[cl]=rep[0].split(' ')[1].replace('...','').replace('>','')
    for cl in cluster_table:
        cluster_table[cl].append(cluster_reps[cl])
    #saving as binary object and cluster metadata table
    with open(os.path.join(cluster_p,f"{cluster_pkl_name}.json"),'w') as save:
        json.dump(cluster_table,save)
    with open(os.path.join(cluster_p,cluster_meta_name),'w') as file:
        file.write(f'Cluster\tRepresentative\tGenotypes\tSegments\tHosts\tCountries\tCollection_date\n')
        for elem in cluster_table:
            file.write(
                f'{elem}\t{cluster_table[elem][5]}\t{cluster_table[elem][0]}\t{cluster_table[elem][1]}\t'
                f'{cluster_table[elem][2]}\t{cluster_table[elem][3]}\t{cluster_table[elem][4]}\n'
            )

if __name__ == '__main__':
    config_file = str(Path(sys.argv[1]).expanduser().resolve())
    project_root = Path(__file__).resolve().parent

    #check if config is available
    if not os.path.exists(config_file):
        print(f'Error: {config_file} not found.')
        sys.exit(1)

    #read config file
    config = configparser.ConfigParser()
    config.read(config_file)
    if 'update_reports' not in config['Paths']:
        config['Paths']['update_reports'] = 'update/reports'
    #check if required paths exist
    for path in ['samples', 'runs', 'references', 'reports', 'update_reports', 'logs', 'blast_database', 'cluster_database', 'metadata','update']:
        full_path = _resolve_config_path(config['Paths'][path], project_root)
        os.makedirs(full_path, exist_ok=True)
        # Update config object in memory so subsequent calls use the full path
        config['Paths'][path] = full_path
    #if os.environ.get('CONDA_PREFIX') != config['Functions']['conda_env']:
        #print('Error: conda environment not activated.')
        #sys.exit(1)
    #check metadata and references files
    metadata=config['Filenames']['metadata']
    references=config['Filenames']['ref_db']
    geo_json=config['Filenames']['geo']
    taxa_json=config['Filenames']['taxa']
    metadata_source = project_root / metadata
    references_source = project_root / references
    if not metadata_source.exists():
        print(f'Error: {metadata} not found.')
        sys.exit(1)
    else:
        shutil.move(str(metadata_source), os.path.join(config['Paths']['metadata'], metadata))
    for filename in [geo_json, taxa_json]:
        metadata_target = os.path.join(config['Paths']['metadata'], filename)
        source_path = project_root / filename
        if source_path.exists():
            shutil.move(str(source_path), metadata_target)
            continue
        if os.path.exists(metadata_target):
            print(f'{filename} already exists in metadata. Skipping move.')
            continue
        print(f'Error: {filename} not found.')
        sys.exit(1)
    if not references_source.exists():
        with open(references_source,'a') as f:
            os.utime(references_source, None)
            shutil.move(str(references_source), os.path.join(config['Paths']['references'], references))
    else:
        shutil.move(str(references_source), os.path.join(config['Paths']['references'], references))
    print('Files moved successfully.')
    print('Creating blast database...')
    #create blast database
    if not os.path.exists(os.path.join(config["Paths"]["blast_database"],f'{config["Filenames"]["l_blast"]}.ndb')):
        subprocess.run(['makeblastdb', '-in', f'{config["Filenames"]["l_blast"]}.fasta', '-dbtype', 'nucl', '-out', os.path.join(config["Paths"]["blast_database"],config["Filenames"]["l_blast"])])
    #creating cluster database
    print('Creating cluster database...')
    if not os.path.exists(os.path.join(config["Paths"]["cluster_database"],config["Filenames"]["cluster_clstr"])):    
        subprocess.run(['cd-hit-est', '-i', f'{config["Filenames"]["l_blast"]}.fasta', '-o', os.path.join(config["Paths"]["cluster_database"],config["Filenames"]["cluster"]), '-c', config["CD-HIT"]["identity"], '-g', '1', '-M', config["CD-HIT"]["memory"]])
    #creating cluster metadata file
    print('Creating cluster metadata...')
    if not os.path.exists(os.path.join(config["Paths"]["cluster_database"],config["Filenames"]["cluster_metadata"])):
        cluster_charact(metadata,config['Paths']['metadata'],config['Filenames']['cluster_clstr'],\
            config['Paths']['cluster_database'],config['Filenames']['cluster_pkl'].replace('.json','').replace('.pkl',''),config['Filenames']['cluster_metadata'])
    #creating representative blast database
    print('Creating representative blast database...')
    if not os.path.exists(os.path.join(config["Paths"]["blast_database"], config["Filenames"]["cluster"])):    
        subprocess.run(['makeblastdb', '-in', os.path.join(config["Paths"]["cluster_database"],config["Filenames"]["cluster"]), '-dbtype', 'nucl', '-out', os.path.join(config["Paths"]["blast_database"], config["Filenames"]["cluster"])])
