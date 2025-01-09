from gb_utils import fetch_genbank_list, append_genbank_from_list, get_gb, create_lookup, write_gb,gb_to_fasta
import pandas as pd
import sys
import os

if __name__=='__main__':
    email=sys.argv[3]
    report_path=f'reports/{sys.argv[1]}'
    db_path=f'references/{sys.argv[2]}'

update_local_db=False

#opening report path and fetching references
report=pd.read_table(report_path,index_col=False)
references=report['REPRESENTATIVE'].to_list()
references=set(references)

found=[]
missing=[]

#checking if local_db has entries
records=create_lookup(db_path)
for ref in references:
    if ref in records:
        found.append(ref)
    else:
        missing.append(ref)

#getting missing records
remote=fetch_genbank_list(missing,email)
#getting local records
local=get_gb(db_path, found)
#creating output files
output=remote+local
#writing output files (.gb and .fasta)
outname=(report_path.split('/')[1].split('.')[0])+'_references'
write_gb(output,f'references/{outname}.gb')
gb_to_fasta(f'references/{outname}.gb',f'references/{outname}.fasta')

if update_local_db:
    append_genbank_from_list(db_path,remote)
