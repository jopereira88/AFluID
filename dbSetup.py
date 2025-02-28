import pandas as pd
from gb_utils import fetch_fasta, filter_fasta_by_accession
from flu_utils import pkl_save,get_ncbi_accessions,pkl_load
import sys
import os

if __name__=='__main__':
    metadata_path=sys.argv[1]
    output_name=sys.argv[2]

if os.path.exists(metadata_path):
    end=False
    # start by opening metadata and retrieving sequence accession numbers
    metadata=pd.read_csv(metadata_path,index_col=False)
    accessions=metadata['ACCESSION'].to_list()

    while not end:
        print('Sequence data downloader\nSelect option:')
        print('1. I have a NCBI sequence multifasta.')
        print("2. I don't have a NCBI Sequence multifasta.")
        print('3. I have a missed accession .pkl file and want to attempt downloading')
        print('4. Exit')
        opt=int(input('Enter option:  '))
        if opt in (1,2,3,4):
            if opt == 1:
                fasta_path=input('Enter fasta filepath: ')
                if os.path.exists(fasta_path):
                    owned=get_ncbi_accessions(fasta_path)
                    db=set(accessions)
                    missing=db-owned
                    if not missing:
                        filter_fasta_by_accession(fasta_path,accessions,output_name)
                    else:
                        missed=list(missed)
                        print(f'The following accessions are missing from the fasta file: {missed}')
                    end=True
                else:
                    print('Invalid path') 
            elif opt == 2:
                # getting the fasta file
                print('WARNING: Downloading a multi-fasta file from NCBI can take a long time\
                      and may fail due to network issues. If you have a fasta file, select option 1.')
                print('Do you want to download? ')
                choice=input('(y/n): ')
                if choice.lower()=='y':
                    email=input('Enter NCBI ENTREZ email: ')
                    missed=fetch_fasta(accessions,f'{output_name}.fasta',email=email)
                    if missed:
                        pkl_save(missed,'failed_accessions')
                    end=True
                else:
                    opt=4
            elif opt==3:
                print('WARNING: Downloading a multi-fasta file from NCBI can take a long time\
                      and may fail due to network issues. If you have a fasta file, select option 1.')
                print('Do you want to download? ')
                choice=input('(y/n): ')
                if choice.lower()=='y':
                    missed=pkl_load('failed_accessions.pkl')
                    email=input('Enter NCBI ENTREZ email: ')
                    fetch_fasta(missed,f'{output_name}.fasta',email=email)
                else:
                    opt=4
            else:
                print('Exiting program')
                end=True
        else:
            print('Invalid option\nSelect a valid option')

else:
    print('Metadata file not found, exiting program')


