#! /bin/python3

from Bio import SeqIO, Entrez
from flu_utils import headers_from_mult_fas, pkl_save, pkl_load
from time import sleep

Entrez.email =  'jlgp1988@gmail.com'

def fetch_genbank(access, file_name, db='nucleotide',start_acc=None,write_mode='a'):
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
        with open(file_name, write_mode) as output:
            if start_acc:
                for i, acc in enumerate(access[access.index(start_acc):]):
                    if i//5000 == 0:
                        sleep(300)
                        with Entrez.efetch(db=db,id=acc, rettype='gb',retmode='text') as handle:
                            data=handle.read()
                            output.write(data)
                    else:
                        with Entrez.efetch(db=db,id=acc, rettype='gb',retmode='text') as handle:
                            data=handle.read()
                            output.write(data)    
            else:
                for i, acc in enumerate(access):
                    if i//5000==0:
                        sleep(300)
                        with Entrez.efetch(db=db,id=acc, rettype='gb',retmode='text') as handle:
                            data=handle.read()
                            output.write(data)
                    else:
                        with Entrez.efetch(db=db,id=acc, rettype='gb',retmode='text') as handle:
                            data=handle.read()
                            output.write(data)
    
    except Exception as e:
        print(f'Error fetching data for accession numbers: {access}, Error: {e}')
        failed.append(access)
    
    return failed

def filter_genbank_by_accession(input_file, output_file, accession_list, verbose=False):
    """
    Filters a multi-entry GenBank file to keep only records with accessions in a specified list.
    Also reports accessions in the list that are not found in the GenBank file.
    
    Args:
        input_file (str): Path to the input multi-entry GenBank file.
        output_file (str): Path to the output GenBank file with filtered entries.
        accession_list (list of str): List of accession numbers (with versions) to keep.
    
    Returns:
        list of str: Accessions from the list that were not found in the GenBank file.
    """
    # Convert accession list to a set for fast lookup
    accession_set = set(accession_list)
    found_accessions = set() 

    # Parse through the GenBank file and filter entries
    matching_records = []
    for record in SeqIO.parse(input_file, "genbank"):
        # Construct the versioned accession (e.g., 'NC_000852.1')
        gb_accession = f"{record.id}"

        # Check if the current record's versioned accession is in the list
        if gb_accession in accession_set:
            matching_records.append(record)
            found_accessions.add(gb_accession)  # Track as found
    # Determine which accessions were not found in the file
    
    missing_accessions = accession_set - found_accessions

    # Write the filtered records to the output file
    
    SeqIO.write(matching_records, output_file, "genbank")
    if verbose: 
        print(f"Filtered GenBank file saved to {output_file} with {len(matching_records)} records.")

    # Output the missing accessions
    if verbose:    
        if missing_accessions:
            print(f"Accessions not found in the file: {', '.join(missing_accessions)}")
        else:
            print("All accessions from the list were found in the file.")

    return list(missing_accessions)

def append_genbank_from_gb(main_file, other_file, verbose=False):
    '''
    
    '''
    #Create a set with the two gb record ID's
    main_records=set()
    other_records=set()
    if verbose:
        print(f'Parsing {main_file}')
    for record in SeqIO.parse(main_file,'genbank'):
        rec=f'{record.id}'
        main_records.add(rec)
    
    if verbose:
        print(f'Parsing {other_file}')
    for record in SeqIO.parse(other_file,'genbank'):
        rec=f'{record.id}'
        other_records.add(rec)
    
    # evaluating for number of missing records
    if verbose:
        print('Evaluating for missing records')
    
    infile=list(SeqIO.parse(main_file,'genbank'))
    missing_records=other_records-main_records
    
    if list(missing_records)==[]:
        print(f'All records from {other_file} in {main_file}')
    else:
        if verbose:
            print(f'Appending missing records to {main_file}')
        outfile=[]
        for record in SeqIO.parse(other_file,'genbank'):
            if record.id in missing_records:
                outfile.append(record)
        
        all_records= infile + outfile
        SeqIO.write(all_records, main_file,'genbank')
        if verbose:
            print(f'Appended {len(outfile)} new records to {main_file}')

def check_missing(gb_file,accession_list):
    '''
    
    '''
    # Mines genbank ID's into a set
    accession_set=set(accession_list)
    gb_set=set()
    for record in SeqIO.parse(gb_file,'genbank'):
        rec=f'{record.id}'
        gb_set.add(rec)
    #missing id's
    missing=accession_set-gb_set
    if list(missing)==[]:
        print('No missing records')
    else:
        return list(missing)


accessions=headers_from_mult_fas(['sequencesDnaInf.fasta'],only_name=True)
missing=pkl_load('acc_to_entrez.pkl')
#filter_genbank_by_accession('sequence.gb','sequencesDnaInf_test.gb',accessions)
failed=fetch_genbank(missing,'sequencesDnaInf_1.gb')
#append_genbank_from_gb('sequencesDnaInf_test.gb','sequencesDnaInf_1.gb',verbose=True)

pkl_save(failed,'acc_to_entrez')