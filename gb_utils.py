#! /bin/python3

from Bio import SeqIO, Entrez
from io import StringIO
from time import sleep



def fetch_genbank(access, file_name, email, db='nucleotide',start_acc=None,write_mode='a'):
    '''
    Fetches GenBank files from NCBI given a list of accession numbers and saves them into a file
    Args:
    access: list of accession numbers
    file_name: str, name of the output file
    db: str, database to fetch from (default is 'nucleotide')
    email: str, Entrez email
    Returns: list of failed files
    '''
    Entrez.email = email

    failed=[]
    try:
        with open(file_name, write_mode) as output:
            if start_acc:
                for i, acc in enumerate(access[access.index(start_acc):]):
                    if i%5000 == 0 and i!=0:
                        sleep(60)
                        with Entrez.efetch(db=db,id=acc, rettype='gb',retmode='text') as handle:
                            data=handle.read()
                            output.write(data)
                    else:
                        with Entrez.efetch(db=db,id=acc, rettype='gb',retmode='text') as handle:
                            data=handle.read()
                            output.write(data)    
            else:
                for i, acc in enumerate(access):
                    if i%5000 == 0 and i!=0:
                        sleep(60)
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
    if failed != []:
        print('Accessions that were not retrieved:')
        for elem in failed:
            print(elem)

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
def gb_to_fasta(input_file,output_file):
    '''
    Converts a gb file to fasta file (header+full sequence)

    '''
    with open(output_file, "w") as fasta_file:
        SeqIO.write(SeqIO.parse(input_file, "genbank"), fasta_file, "fasta")    

def gb_feature_to_fasta(input_file, output_file, feature_type='CDS'):
    '''
    Converts a gb file into a fasta file
    param input_file(str): path for input gb file
    param output_file(str): path and name of the output fasta file
    param feture_type(str): type of feature to mine the sequence CDS/gene
    (Default = 'CDS')
    '''
    with open(output_file, "w") as fasta_file:
        for record in SeqIO.parse(input_file, "genbank"):
            record_id=f'{record.id}'
            for feature in record.features:
                if feature.type == feature_type:  # Look for CDS features
                # Extract the sequence
                    feature_seq = feature.location.extract(record.seq)
                # Create a FASTA record
                    fasta_file.write(f">{record_id} | {feature_type}_{feature.qualifiers.get('gene',['unknown'])[0]}\n")
                    fasta_file.write(str(feature_seq) + "\n")

def get_gb(gb_file,accessions):
    '''
    Fetches a genbank record from a gb file by acession numbers

    '''
    outrecords=[]
    for record in SeqIO.parse(gb_file,'genbank'):
        if record.id in accessions:
            outrecords.append(record)
    return outrecords

def create_lookup(gbfile):
    ids=set()
    for record in SeqIO.parse(gbfile,'genbank'):
        ids.add(record.id)
    return ids

def is_record(gbfile,accession):
    '''
    Parses a gb file to check if an accession exists as a record in a gb file
    '''
    for record in SeqIO.parse(gbfile,'genbank'):
        if record.id == accession:
            return True
    return False
def fetch_genbank_list(access, email, db='nucleotide'):
    '''
    Fetches GenBank files from NCBI given a list of accession numbers and saves them into a list
    access: list of accession numbers
    db: str, database to fetch from (default is 'nucleotide')
    email: str, Entrez email
    Returns: list of failed files
    '''
    Entrez.email = email 
    failed=[]
    output=[]
    try:
        for i, acc in enumerate(access):
            if i%5000 == 0 and i!=0:
                sleep(60)
                with Entrez.efetch(db=db,id=acc, rettype='gb',retmode='text') as handle:
                    record=handle.read()
                    output.append(record)
            else:
                with Entrez.efetch(db=db,id=acc, rettype='gb',retmode='text') as handle:
                    record=handle.read()
                    output.append(record)
    except Exception as e:
        print(f'Error fetching data for accession numbers: {access}, Error: {e}')
        failed.append(access)
        if failed != []:
            print('Accessions that were not retrieved:')
            for elem in failed:
                print(elem)
    parsed_output=[]
    for record in output:
        record_io=StringIO(record)
        for seq_record in SeqIO.parse(record_io,'genbank'):
            parsed_output.append(seq_record)
    return parsed_output

def write_gb(records,output):
    with open(output,'w') as gbfile:
        SeqIO.write(records,gbfile,'genbank')

def append_genbank_from_list(main_file, records, verbose=False):
    '''
    
    '''
    #Create a set with the two gb record ID's
    main_records=set()
    other_records=set()

    for record in SeqIO.parse(main_file,'genbank'):
        rec=f'{record.id}'
        main_records.add(rec)
    for record in records:
        other_records.add(record.id)
    
    infile=list(SeqIO.parse(main_file,'genbank'))
    missing_records=other_records-main_records
    
    if list(missing_records)==[]:
        print(f'All records from the record list in {main_file}')
    else:
        if verbose:
            print(f'Appending missing records to {main_file}')
        outfile=[]
        for record in records:
            if record.id in missing_records:
                outfile.append(record)
        
        all_records= infile + outfile
        SeqIO.write(all_records, main_file,'genbank')
        if verbose:
            print(f'Appended {len(outfile)} new records to {main_file}')

def fetch_fasta(access, file_name, email, db='nucleotide',start_acc=None,write_mode='a'):
    '''
    Fetches Fasta files from NCBI given a list of accession numbers and saves them into a file
    Args:
    access: list of accession numbers
    file_name: str, name of the output file
    db: str, database to fetch from (default is 'nucleotide')
    email: str, Entrez email
    Returns: list of failed files
    '''
    Entrez.email = email

    failed=[]
    try:
        with open(file_name, write_mode) as output:
            if start_acc:
                for i, acc in enumerate(access[access.index(start_acc):]):
                    if i%5000 == 0 and i!=0:
                        sleep(60)
                        with Entrez.efetch(db=db,id=acc, rettype='fasta',retmode='text') as handle:
                            data=handle.read()
                            output.write(data)
                    else:
                        with Entrez.efetch(db=db,id=acc, rettype='fasta',retmode='text') as handle:
                            data=handle.read()
                            output.write(data)    
            else:
                for i, acc in enumerate(access):
                    if i%5000 == 0 and i!=0:
                        sleep(60)
                        with Entrez.efetch(db=db,id=acc, rettype='fasta',retmode='text') as handle:
                            data=handle.read()
                            output.write(data)
                    else:
                        with Entrez.efetch(db=db,id=acc, rettype='fasta',retmode='text') as handle:
                            data=handle.read()
                            output.write(data)
    
    except Exception as e:
        print(f'Error fetching data for accession numbers: {access}, Error: {e}')
        failed.append(access)
    if failed != []:
        print('Accessions that were not retrieved:')
        for elem in failed:
            print(elem)
        return failed
def filter_fasta_by_accession(infile,access,outfile):
    outputs=[]
    with open(infile,'r') as handle:
        records=SeqIO.parse(handle,'fasta')
        for record in records:
            if record.id in access:
                outputs.append(record)
    with open(outfile,'w') as filtered:
        for elem in outputs:
            SeqIO.write(elem,filtered,'fasta')

def remove_records(fasta_file,accessions,output_file):
    '''
    Removes records from a fasta file given a list of accessions
    '''
    records=[]
    with open(fasta_file,'r') as handle:
        for record in SeqIO.parse(handle,'fasta'):
            if record.id not in accessions:
                records.append(record)
    with open(output_file,'w') as output:
        for record in records:
            SeqIO.write(record,output,'fasta')