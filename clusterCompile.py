#!/usr/bin/python3
'''This script creates a report file obtained after assigment of clustered sequences

Outputs a .txt tabular file with sample information, %ID against the cluster_rep_seq,
cluster, cluster_rep_seq accession number, genotypes, segments and hosts

Used as a part of the bestCluster.sh script'''

from flu_utils import pkl_load
import sys
import warnings
warnings.filterwarnings("ignore", message=".*does not match any known type.*")

# Variable assignment as script args
if __name__ == '__main__':
    s_dict=sys.argv[1]
    cl_assign=sys.argv[2]
db_path='cluster_db/'
reports_path='reports/'
dict_sample=pkl_load(s_dict)
assign=open(cl_assign,'r').readlines()[1:]
dict_clust=pkl_load(f'{db_path}dict_clusters.pkl')

# Creates sample_assign dict {sample:%ID}
sample_assign={}
for line in assign:
    line=line.strip()
    line=line.split(';')
    if len(line)>2:
        sample_assign[line[1]]=line[2]
    else:
        sample_assign[line[1]]='NA'

#Creates and populates the compiler dict
#{sample:[%ID,Cluster,genotypes,segment,hosts,Cluster_rep]}
compiler={}
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

#Uses compiler to output the report file
output_name=cl_assign.split('/')[2]
output_name=output_name.replace('.assign','')

with open(f'{reports_path}{output_name}_clust.txt','w') as report:
    report.write('Sample_access\t%ID\tCluster\tCluster_rep\tGenotypes\tSegment\tHosts\n')
    for sample in compiler:
        if len(compiler[sample])>2:
           report.write(f'{sample}\t{compiler[sample][0]}\t\
                        {compiler[sample][1]}\t{compiler[sample][5]}\t\
                        {compiler[sample][2]}\t{compiler[sample][3]}\t{compiler[sample][4]}\n')
        elif len(compiler[sample])==0:
            continue  
        else:
            report.write(f'{sample}\t{compiler[sample][0]}\t{compiler[sample][1]}\n')

print(f"Compiled report generated in {reports_path}{output_name}_clust.txt")