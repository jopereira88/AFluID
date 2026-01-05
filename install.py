import os
import sys
import configparser
import shutil
import subprocess
from flu_utils import metadata_dict,parse_clstr
from collections import Counter
import pickle

def cluster_charact(metadata,metadata_p,clstr,cluster_p,cluster_pkl_name,cluster_meta_name):
    meta=metadata_dict(os.path.join(metadata_p,metadata),[9,10,12,11],sep=';')
    clusters=parse_clstr(os.path.join(cluster_p,clstr))
    #getting metadata annotations per cluster
    cluster_table={}
    for key in clusters:
        cluster_table[key]=[Counter(),Counter(),Counter(),Counter()] #genotype,segment,host,country
        for value in clusters[key]:
            cluster_table[key][0][meta[value][0]]+=1
            cluster_table[key][1][meta[value][1]]+=1
            cluster_table[key][2][meta[value][2]]+=1
            cluster_table[key][3][meta[value][3]]+=1
    for elem in cluster_table:
        cluster_table[elem][0]=dict(cluster_table[elem][0])
        cluster_table[elem][1]=dict(cluster_table[elem][1])
        cluster_table[elem][2]=dict(cluster_table[elem][2])
        cluster_table[elem][3]=dict(cluster_table[elem][3])
    #getting cluster representatives
    clusters=parse_clstr(os.path.join(cluster_p,clstr),access_only=False)
    cluster_reps={}
    for cl in clusters:
        rep=[item for item in clusters[cl] if '*' in item]
        cluster_reps[cl]=rep[0].split(' ')[1].replace('...','').replace('>','')
    for cl in cluster_table:
        cluster_table[cl].append(cluster_reps[cl])
    #saving as binary object and cluster metadata table
    with open(os.path.join(cluster_p,cluster_pkl_name),'wb') as save:
        pickle.dump(cluster_table,save)
    with open(os.path.join(cluster_p,cluster_meta_name),'w') as file:
        file.write(f'Cluster\tRepresentative\tGenotypes\tSegments\tHosts\tCountries\n')
        for elem in cluster_table:
            file.write(f'{elem}\t{cluster_table[elem][4]}\t{cluster_table[elem][0]}\t{cluster_table[elem][1]}\t{cluster_table[elem][2]}\t{cluster_table[elem][3]}\n')

if __name__ == '__main__':
    config_file=sys.argv[1]

    #check if config is available
    if not os.path.exists(config_file):
        print(f'Error: {config_file} not found.')
        sys.exit(1)

    #read config file
    config = configparser.ConfigParser()
    config.read(config_file)
    cwd=os.getcwd()
    #check if required paths exist
    for path in ['samples', 'runs', 'references', 'reports', 'logs', 'blast_database', 'cluster_database', 'metadata']:
        full_path = os.path.abspath(os.path.join(cwd, config['Paths'][path]))
        os.makedirs(full_path, exist_ok=True)
        # Update config object in memory so subsequent calls use the full path
        config['Paths'][path] = full_path
    #if os.environ.get('CONDA_PREFIX') != config['Functions']['conda_env']:
        #print('Error: conda environment not activated.')
        #sys.exit(1)
    #check metadata and references files
    metadata=config['Filenames']['metadata']
    references=config['Filenames']['ref_db']
    if not os.path.exists(metadata):
        print(f'Error: {metadata} not found.')
        sys.exit(1)
    else:
        shutil.move(metadata,os.path.join(config['Paths']['metadata'],metadata))
    if not os.path.exists(references):
        with open(references,'a') as f:
            os.utime(references,None)
            shutil.move(references,os.path.join(config['Paths']['references'],references))
    else:
        shutil.move(references,os.path.join(config['Paths']['references'],references))
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
            config['Paths']['cluster_database'],config['Filenames']['cluster_pkl'].replace('.pkl',''),config['Filenames']['cluster_metadata'])
    #creating representative blast database
    print('Creating representative blast database...')
    if not os.path.exists(os.path.join(config["Paths"]["blast_database"], config["Filenames"]["cluster"])):    
        subprocess.run(['makeblastdb', '-in', os.path.join(config["Paths"]["cluster_database"],config["Filenames"]["cluster"]), '-dbtype', 'nucl', '-out', os.path.join(config["Paths"]["blast_database"], config["Filenames"]["cluster"])])
