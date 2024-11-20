'''
This script conforms a fasta file for flumut analysis
'''
#!/usr/bin/python3
from flu_utils import pkl_load,seq_get,int_to_iupac, dict_to_fasta
import sys


#fasta file
if __name__ == '__main__':
    fasta=sys.argv[1]

# Paths
sample_path='samples/'
run_path='runs/'

#load ID report
report=pkl_load(f'{run_path}to_flumut.pkl')
#creating list of sequences to filter
#create fasta file dict
flt_fasta=seq_get(f'{sample_path}{fasta}')

#mining segments
for key in report:
    report[key][3]=eval(str(report[key][3]))
    report[key][3]=[f'>{key}|{int_to_iupac(int(i))}' for i in report[key][3] if len(report[key][3])==1]

out_fasta={}
for key in flt_fasta:
    try:    
        out_fasta[report[key.replace('>','')][3][0]]=flt_fasta[key]
    except KeyError:
        continue
dict_to_fasta(out_fasta, f'{sample_path}to_flumut')