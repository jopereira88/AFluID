#!/usr/bin/python3
from flu_utils import pkl_load, pkl_save
import sys
import os
from metadata_utils import BlastReportTable, BlastClustReportTable



###################################################################
runs_dir='runs/'
samples_dir='samples/'
to_report=pkl_load(f'{runs_dir}assigned.pkl')
if __name__=='__main__':
    filename1=sys.argv[1]
    filename2=sys.argv[2]
    outname=sys.argv[3] 

#filename1='to_blast_bclust_report.txt'
#filename2='to_reblast_report.txt'
#outname='test_report.txt'
report_path='reports/'
to_remote=[]
if os.path.exists(f'{report_path}{filename1}'):
    bclust_report=BlastClustReportTable(f'{report_path}{filename1}')
    bclust_report.convert_to_prop('genotype')
    bclust_report.convert_to_prop('host')
    for key in bclust_report.data:
        if key not in to_report and bclust_report.data[key][2]!='NA':
            to_report[key]=[]
            to_report[key].append(bclust_report.data[key][0])
            to_report[key].append(bclust_report.data[key][1])
            to_report[key].append(bclust_report.data[key][2])
            to_report[key].append(bclust_report.data[key][3])
            to_report[key].append(bclust_report.data[key][4])
            to_report[key].append(bclust_report.data[key][5])
            to_report[key].append('C-BLAST')
    if os.path.exists(f'{report_path}{filename2}'):
        bblast_report=BlastReportTable(f'{report_path}{filename2}')
        for key in bblast_report.data:
            print(bblast_report.data[key])
            if bblast_report.data[key][1]=='NA':
                to_remote.append(key)
            elif key not in to_report and bblast_report.data[key][1]!='NA':
                to_report[key]=[]
                to_report[key].append(bblast_report.data[key][0])
                to_report[key].append('UNCLUSTERED')
                to_report[key].append(bblast_report.data[key][1])
                to_report[key].append(bblast_report.data[key][2].replace(' ',''))
                to_report[key].append(bblast_report.data[key][3].replace(' ',''))
                to_report[key].append(bblast_report.data[key][4])
                to_report[key].append('L-BLAST')
mappings=pkl_load(f'{samples_dir}{outname}_mappings.pkl')

pkl_save(to_report,f'{runs_dir}to_flumut')

with open(f'{report_path}{outname}_ID_report.txt','w') as report:
        report.write('SAMPLE_NAME\tREPRESENTATIVE\tCLUSTER\t%ID\tSEGMENT\tGENOTYPE\tHOST\tASSIGNED_BY\n')
        for key in to_report:
            mapped=mappings[f'>{key}']
            report.write(f'{key}\t{to_report[key][0]}\t\
                         {to_report[key][1]}\t{to_report[key][2]}\t{to_report[key][3]}\t\
                            {to_report[key][4]}\t{to_report[key][5]}\t{to_report[key][6]}\n')
if len(to_remote)>0:
    print(f'Following samples were unassigned by local BLAST: {len(to_remote)}\nSample accession numbers:')
    for item in to_remote:
        print(item)