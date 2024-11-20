'''
Assigns sample names to Seq values
'''
from flu_utils import pkl_load
import sys
import pandas as pd

samples_dir='samples/'
reports_dir='reports/'
if __name__ == '__main__':
    filename=sys.argv[1]

mappings=pkl_load(f'{samples_dir}{filename}_mappings.pkl')
mapping_dict={key.replace('>',''):value for key , value in mappings.items()}
mk_report=pd.read_table(f'{reports_dir}{filename}_markers.tsv')
mk_report=mk_report.replace(mapping_dict)
mut_report=pd.read_table(f'{reports_dir}{filename}_mutations.tsv')
mut_report=mut_report.replace(mapping_dict)
mk_report.to_csv(f'{reports_dir}{filename}_markers.tsv', sep='\t', index=False)
mut_report.to_csv(f'{reports_dir}{filename}_mutations.tsv', sep='\t', index=False)