#!/usr/bin/python3
from flu_utils import sp_blastn, headers_from_mult_fas, seq_filter_get, dict_to_fasta
from metadata_utils import ClusterMetadata
import sys
#import time
#start_time=time.time()

#Paths
metadata_path='cluster_db/cluster_desc.txt'
data_path='samples/'
db_path='blast_db/infDNAClusters'
run_path='runs/'
output_path='reports/'

if __name__ == '__main__':
    fasta_file=sys.argv[1]

#fasta_file='format_EISN_INF_EQA_coded.fasta' #for direct script access

#priming sample names
access=headers_from_mult_fas([f'{data_path}{fasta_file}'],only_name=True)
queries={key:[] for key in access}
print('Starting BLASTn')
sp_blastn(f'{data_path}{fasta_file}',db_path,f'{run_path}{fasta_file.split(".")[0]}_run',createdb=False,silent=True,maxtargetseqs=7)
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

metadata=ClusterMetadata(metadata_path)
report={}
to_reblast=set()
for i in queries:
        report[i]={}
        if len(queries[i])>0:
            for j in queries[i]:
                if float(j[1])<90.0:
                    to_reblast.add(i)
                else:
                    report[i][j[0]]=[j[1]]
        else:
            report[i]['Unassigned']='NA'
            to_reblast.add(i)
for item in to_reblast:
    if item in report.keys():
        to_reblast.discard(item)

#removing clusters 27 and 22681 from the report
clusters_reps={}
for cluster in metadata.data:
    for i in report.values():
        for key in i.keys():
            if metadata.data[cluster][metadata.headers['REPRESENTATIVE']-1]==key:
                clusters_reps[key]=cluster
keys_to_remove=set()
for key in clusters_reps:
    if clusters_reps[key] == '>Cluster 27' or clusters_reps[key] == '>Cluster 22681':
        keys_to_remove.add(key)
for key in keys_to_remove:
    for i in report:
        if report[i].get(key):
            del report[i][key]
for key in report:
    for rep in clusters_reps:
        if rep in report[key]:
            report[key][rep].append(clusters_reps[rep])
            report[key][rep].append(metadata.data[clusters_reps[rep]][metadata.headers['SEGMENTS']-1])
            report[key][rep].append(metadata.data[clusters_reps[rep]][metadata.headers['GENOTYPES']-1])
            report[key][rep].append(metadata.data[clusters_reps[rep]][metadata.headers['HOSTS']-1])          

#outputing report
with open(f'{output_path}{fasta_file.split(".")[0]}_bclust_report.txt','w') as repo:
    repo.write(f'Sample ID\tCluster_Rep\tCluster\t%identity\tSegment\tGenotype\tHost\n')
    for item in report:
        for j in report[item]:
            if report[item][j]!='NA':
                repo.write(f'{item}\t{j}\t{report[item][j][1]}\t{report[item][j][0]}\t\
                        {report[item][j][2]}\t{report[item][j][3]}\t{report[item][j][4]}\n')
            else:
                repo.write(f'{item}\t{j}\t{report[item][j]}\n')

print(f'Outputted report to path:{output_path}{fasta_file.split(".")[0]}_bclust_report.txt')
if list(to_reblast) != []:
   to_reblast_dict=seq_filter_get(f'{data_path}{fasta_file}',to_reblast)
   dict_to_fasta(to_reblast_dict,f'{data_path}to_reblast') 

#end_time=time.time()
#exec_time=end_time-start_time
#print(f'Completed in {exec_time} seconds')