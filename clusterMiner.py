#!/usr/bin/python3

#Importing modules
from metadata_utils import ClusterReportTable
from flu_utils import pkl_save, seq_filter_get, dict_to_fasta
import sys
import warnings
warnings.filterwarnings("ignore", message=".*does not match any known type.*")

#Parsing arguments
if __name__=='__main__':
    filename=sys.argv[1]
    sample_fasta=sys.argv[2]


#setting paths
report_path='reports/'
fasta_path='samples/'
runs_dir='runs/'

#Initiating a cluster report table object
clust_rep=ClusterReportTable(f'{report_path}{filename}')

#Converting counts into percentages
clust_rep.convert_to_prop('genotypes')
clust_rep.convert_to_prop('hosts')

#Getting information on hemagglutinin and neuraminidase genotypes
hs=clust_rep.extract_h()
ns=clust_rep.extract_n()

#Filtering data to blast
to_blast=[]
for key in clust_rep.data:
    if clust_rep.data[key][0]=='NA': #removes unassigned sequences
        to_blast.append(key)
    elif len(clust_rep.data[key][4]) > 1: #removes sequences from segment-ambiguous clusters
        to_blast.append(key)
    #removes sequences that have genotype ambiguous clusters
    elif clust_rep.data[key][4] == '4': 
        if len(hs[key]) > 1:
            to_blast.append(key)
    elif clust_rep.data[key][4] =='6':
        if len(ns[key]) > 1:
            to_blast.append(key)
print(to_blast)
#Writing report dict
to_report={}
for key in clust_rep.data:
    if key not in to_blast:
        to_report[key]=[]
        to_report[key].append(clust_rep.data[key][2].replace(' ',''))
        to_report[key].append(
            clust_rep.data[key][1].replace('                        >','>'))
        to_report[key].append(clust_rep.data[key][0].replace('%',''))
        to_report[key].append(clust_rep.data[key][4])
        to_report[key].append(clust_rep.data[key][3])
        to_report[key].append(clust_rep.data[key][5])
        to_report[key].append('CD-HIT')

pkl_save(to_report,f'{runs_dir}assigned')

#Outputting to blast fasta file
to_blast_dict=seq_filter_get(f'{fasta_path}{sample_fasta}',to_blast)
dict_to_fasta(to_blast_dict,f'{fasta_path}to_blast')
