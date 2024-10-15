#!/usr/bin/python3
from flu_utils import pkl_load
import sys
from metadata_utils import BlastReportTable
import warnings
warnings.filterwarnings("ignore", message=".*does not match any known type.*")


###################################################################
runs_dir='runs/'
to_report=pkl_load(f'{runs_dir}assigned.pkl')
if __name__=='__main__':
    filename=sys.argv[1]
    outname=sys.argv[2]
report_path='reports/'
bblast_report=BlastReportTable(f'{report_path}{filename}')
blast_segs=bblast_report.get_segment()
blast_gen=bblast_report.get_genotype()
blast_host=bblast_report.get_host()
to_remote=[]
for key in bblast_report.data:
    if bblast_report.data[key][0]=='NA':
        to_remote.append(key)
    else:
        to_report[key]=[]
        to_report[key].append(bblast_report.data[key][0])
        to_report[key].append(bblast_report.data[key][1])
        to_report[key].append(blast_gen[key])
        to_report[key].append(blast_segs[key])
        to_report[key].append(blast_host[key])
        to_report[key].append('BLAST')

with open(f'reports/{outname}_final_report.txt','w') as report:
    report.write('SAMPLE_NAME\t%ID\tREPRESENTATIVE\tGENOTYPE\tSEGMENT\tHOST\tASSIGNED_BY\n')
    for key in to_report:
        report.write(f'{key}\t{to_report[key][0]}\t\
                     {to_report[key][1]}\t{to_report[key][2]}\t\
                        {to_report[key][3]}\t{to_report[key][4]}\t{to_report[key][5]}\n')