#!/usr/bin/python3
'''This script creates a file with the cluster assignment after 
the cd-hit-est-2d run 

Outputs a .csv file with cluster, representative sequence, genotypes, 
segments and hosts

Used as a part of the bestCluster.sh script'''

from flu_utils import parse_clstr,pkl_save,headers_from_mult_fas
import sys
import warnings
warnings.filterwarnings("ignore", message=".*does not match any known type.*")


#Initialising script arguments
if __name__ == '__main__':
    report=sys.argv[1]
    headers_file=sys.argv[2]

#report='runs/rand_100_0.clstr'
#headers_file='samples/rand_100_0.fasta'

#Parsing the .clstr file and creating the assignment dict
clusters=parse_clstr(report,access_only=False)
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
                assigned.append(name)
                report[cl].append((name,id))

#Parsing .fasta file to find unassigned samples
headers=headers_from_mult_fas([headers_file],only_name=True)
unassigned=[item.strip() for item in headers if item.strip() not in assigned]
report['Unassigned']=unassigned
#Creating a samples dict to allow future parsing of sample/cluster relation
sample_dict={}
for header in headers:
    header=header.replace('|','')
    header=header.replace(' ','')
    header=header.replace('>','')  
    sample_dict[header]=''
for key in report:   
    for i in range(len(report[key])):
        if type(report[key][i])==tuple:
            sample_dict[report[key][i][0]]=key
        else:
            sample_dict[report[key][i]]=f'>{key}'

#Creating a .pkl file to export the dict to a future script
pkl_save(sample_dict,'runs/dict_sample')

#Writing the assignment report
output_name=headers_file.split('/')[-1]
output_name=output_name.replace('.fasta','.assign')

with open(f'runs/{output_name}','w') as repo:
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
