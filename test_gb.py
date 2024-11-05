from Bio import SeqIO, Entrez
from flu_utils import headers_from_mult_fas
from progress.bar import Bar

Entrez.email =  'jlgp1988@gmail.com'

def fetch_genbank(access, file_name, db='nucleotide',start_acc=None):
    '''
    Fetches GenBank files from NCBI given a list of accession numbers and saves them into a file
    Args:
    access: list of accession numbers
    file_name: str, name of the output file
    db: str, database to fetch from (default is 'nucleotide')
    Returns: list of failed files
    '''
    failed=[]
    try:
        with Bar('Downloading gbk', fill='#', suffix='%(percent).1f%% - %(eta)ds') as bar:
            with open(file_name, 'w') as output:
                if start_acc:
                    for acc in access[access.index(start_acc):]:
                        with Entrez.efetch(db=db,id=acc, rettype='gb',retmode='text') as handle:
                            data=handle.read()
                            output.write(data)
                            output.write('\n')
                    bar.next()       
                else:
                    for acc in access:
                        with Entrez.efetch(db=db,id=acc, rettype='gb',retmode='text') as handle:
                            data=handle.read()
                            output.write(data)
                            output.write('\n')
                    bar.next()
    
    except Exception as e:
        print(f'Error fetching data for accession numbers: {access}, Error: {e}')
        failed.append(access)
    return failed

accessions=headers_from_mult_fas(['sequencesDnaInf.fasta'],only_name=True)


fetch_genbank(accessions,'sequencesDnaInf_1.gb',start_acc='PP997625.1')