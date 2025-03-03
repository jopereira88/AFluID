import sys
import os
import pandas as pd
from gb_utils import remove_records
from flu_utils import get_ncbi_accessions

if __name__=='__main__':
    metadata_path=sys.argv[1]
    fasta_path=sys.argv[2]
    output_name=sys.argv[3]
    try:
        if os.path.exists(metadata_path):
            # start by opening metadata and retrieving sequence accession numbers
            metadata=pd.read_csv(metadata_path,index_col=False)
            accessions=metadata['ACCESSION'].to_list()
            if os.path.exists(fasta_path):
                owned=get_ncbi_accessions(fasta_path)
                db=set(accessions)
                missing=db-owned
                to_rm=owned-db
                if not missing:
                    remove_records(fasta_path,to_rm,output_name)
                elif missing and to_rm:
                    remove_records(fasta_path,to_rm,output_name)
                    missing=list(missing)
                    print(f'The following accessions are missing from the fasta file: {missing}')
                    with open('missing_accessions.txt','w') as f:
                        for id in missing:
                            f.write(f'{id}\n')
                else:
                    print(f'The following accessions are missing from the fasta file: {missing}')
                    with open('missing_accessions.txt','w') as f:
                        for id in missing:
                            f.write(f'{id}\n')
    except IndexError:
        print("Invalid number of arguments")
        print('Usage:')
        print('python3 dbSetup.py [PATH_TO_METADATA_FILE] [PATH_TO_FASTA_FILE] [OUTPUT_FILE_NAME] ')

