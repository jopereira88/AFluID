#!/usr/bin/python3
from flu_utils import pkl_load
import sys
import os
from metadata_utils import BlastReportTable, BlastClustReportTable



###################################################################
runs_dir='runs/'
to_report=pkl_load(f'{runs_dir}assigned.pkl')
if __name__=='__main__':
    filename1=sys.argv[1]
    filename2=sys.argv[2]
    outname=sys.argv[3]
#filename1='format_EISN_INF_EQA_coded_bclust_report.txt'
#filename2='format_EISN_INF_EQA_coded_blast_report.txt'
#outname='test_report.txt'
report_path='reports/'
to_remote=[]
if os.path.exists(f'{report_path}{filename1}'):
    print('report found 1')
    bclust_report=BlastClustReportTable(f'{report_path}{filename1}')
    bclust_report.convert_to_prop('genotype')
    bclust_report.convert_to_prop('host')
    for key in bclust_report.data:
        if key not in to_report:
            to_report[key]=[]
            to_report[key].append(bclust_report.data[key][0])
            to_report[key].append(bclust_report.data[key][1])
            to_report[key].append(bclust_report.data[key][2])
            to_report[key].append(bclust_report.data[key][3])
            to_report[key].append(bclust_report.data[key][4])
            to_report[key].append(bclust_report.data[key][5])
            to_report[key].append('C-BLAST')
    if os.path.exists(f'{report_path}{filename2}'):
        print('report found 2')
        bblast_report=BlastReportTable(f'{report_path}{filename2}')
        blast_segs=bblast_report.get_segment()
        blast_gen=bblast_report.get_genotype()
        blast_host=bblast_report.get_host()
        for key in bblast_report.data:
            if bblast_report.data[key][0]=='NA':
                to_remote.append(key)
            elif key not in to_report:
                to_report[key]=[]
                to_report[key].append(bblast_report.data[key][0])
                to_report[key].append(bblast_report.data[key][1])
                to_report[key].append('UNCLUSTERED')
                to_report[key].append(bblast_report.data[key][2])
                to_report[key].append(bblast_report.data[key][3])
                to_report[key].append(bblast_report.data[key][4])
                to_report[key].append('L-BLAST')

with open(f'reports/{outname}_final_report.txt','w') as report:
        report.write('SAMPLE_NAME\tREPRESENTATIVE\tCLUSTER\t%ID\tSEGMENT\tGENOTYPE\tHOST\tASSIGNED_BY\n')
        for key in to_report:
            report.write(f'{key}\t{to_report[key][0]}\t\
                         {to_report[key][1]}\t{to_report[key][2]}\t{to_report[key][3]}\t\
                            {to_report[key][4]}\t{to_report[key][5]}\t{to_report[key][6]}\n')
if len(to_remote)>0:
    print(f'Following samples were unassigned by local BLAST: {len(to_remote)}\nSample accession numbers:')
    for item in to_remote:
        print(item)