#!/usr/bin/python3
from flu_utils import sp_blastn, headers_from_mult_fas
from metadata_utils import SequenceMetadata
import sys
#import time
#start_time=time.time()

#Paths
metadata_path='metadata/'
data_path='samples/'
db_path='blast_db/sequencesDnaInf'
run_path='runs/'
output_path='reports/'

if __name__ == '__main__':
    fasta_file=sys.argv[1]
    metadata_file=f'{metadata_path}{sys.argv[2]}'

#fasta_file='to_reblast.fasta' #for direct script access

#priming sample names
access=headers_from_mult_fas([f'{data_path}{fasta_file}'],only_name=True)
queries={key:[] for key in access}
print('Starting BLASTn')
sp_blastn(f'{data_path}{fasta_file}',db_path,f'{run_path}{fasta_file.split(".")[0]}_run',createdb=False,silent=True,maxtargetseqs=7,numthreads=2)
print('Ending BLASTn')
#creating the run dictionary
with open(f'{run_path}{fasta_file.split(".")[0]}_run.txt') as tabular:
    tabular=tabular.readlines()
    for line in tabular:
        line=line.split('\t')
        if line[0] in queries.keys():
            queries[line[0]].append((line[1],line[2]))

#Opening metadata and compiling the report dictionary
#Structure: Sample ID,%identity, bb_access,  bb_org_name, genotype, segment, host

metadata=SequenceMetadata(metadata_file)
report={}
for i in queries:
        report[i]={}
        if len(queries[i])>0:
            for j in queries[i]:
                report[i][j[0].replace(' ','')]=[j[1]]
        else:
            report[i]['Unassigned']='NA'
for key in report:
    for i in report[key]:
        if report[key][i]!='NA':
                report[key][i].append(metadata.data[i][metadata.headers['SEGMENT']-1])
                report[key][i].append(metadata.data[i][metadata.headers['GENOTYPE']-1])
                report[key][i].append(metadata.data[i][metadata.headers['HOST']-1])



#outputing report
with open(f'{output_path}{fasta_file.split(".")[0]}_report.txt','w') as repo:
    repo.write(f'Sample ID\tBb_access\t%identity\tSegment\tGenotype\tHost\n')
    for item in report:
        for j in report[item]:
            if report[item][j]!='NA':
                repo.write(f'{item}\t{j}\t{report[item][j][0]}\t{report[item][j][1]}\t\
                        {report[item][j][2]}\t{report[item][j][3]}\n')
            else:
                repo.write(f'{item}\t{j}\t{report[item][j]}\t{report[item][j]}\t\
                           {report[item][j]}\t{report[item][j]}\n')
           

#end_time=time.time()
#exec_time=end_time-start_time
#print(f'Completed in {exec_time} seconds')