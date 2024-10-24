#!/usr/bin/python3
from flu_utils import sp_blastn, metadata_dict
import sys
import time
import warnings
warnings.filterwarnings("ignore", message=".*does not match any known type.*")


#Paths
metadata_path='metadata/flu_metadata_4.1_flt.csv'
data_path='samples/'
db_path='blast_db/sequencesDnaInf'
run_path='runs/'
output_path='reports/'

if __name__ == '__main__':
    fasta_file=sys.argv[1]

#fasta_file='rand_100_3.fasta' #for direct script access

#priming sample names
queries={}
with open(f'{data_path}{fasta_file}','r') as query:
    for line in query.readlines():
        if '>' in line:
            line=line.split('|')[0]
            line=line.replace('>','')
            line=line.replace(' ','')
            line=line.strip()
            queries[line]=[]
print('Starting BLASTn')
sp_blastn(f'{data_path}{fasta_file}',db_path,f'{run_path}{fasta_file.split(".")[0]}_blastrun',createdb=False,silent=True)
print('Ending BLASTn')
#creating the run dictionary
with open(f'{run_path}{fasta_file.split(".")[0]}_blastrun.txt') as tabular:
    tabular=tabular.readlines()
    for line in tabular:
        line=line.split('\t')
        if line[0] in queries.keys():
            queries[line[0]].append(line[1])
            queries[line[0]].append(line[2])


#Opening metadata and compiling the report dictionary
#Structure: Sample ID,%identity, bb_access,  bb_org_name, genotype, segment, host

metadata=metadata_dict(metadata_path,[1,9,10,12])

report={}
for i in queries:
        report[i]=[]
        if len(queries[i])>0:
            report[i].append(queries[i][0])
            report[i].append(queries[i][1])
        else:
            report[i].append('Unassigned')
            report[i].append('NA')
for key in report:
    if report[key][0]!='Unassigned':
        for value in metadata[report[key][0]]:
            report[key].append(value)

#outputing report
with open(f'{output_path}{fasta_file.split(".")[0]}_report.txt','w') as repo:
    repo.write(f'Sample ID\t%identity\tBb_access\tBb_org_name\tGenotype\tSegment\tHost\n')
    for item in report:
        if report[item][0]!='Unassigned':
            repo.write(f'{item}\t{report[item][1]}\t{report[item][0]}\t{report[item][2]}\t\
                        {report[item][3]}\t{report[item][4]}\t{report[item][5]}\n')
        else:
            repo.write(f'{item}\t{report[item][1]}\t{report[item][0]}\n')
print(f'Outputted report to path:{output_path}{fasta_file.split(".")[0]}_report.txt')
