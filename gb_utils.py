#! /bin/python3

from Bio import SeqIO, Entrez
from io import StringIO
from time import sleep
from typing import Any, Iterable, Optional



def fetch_genbank(access: list[str], file_name: str, email: str, db: str = 'nucleotide',start_acc: Optional[str] = None,write_mode: str = 'a') -> None:
    """
    Fetch GenBank records from NCBI and append them to a file.

    Parameters
    ----------
    access : list
        Accession numbers to fetch.
    file_name : str
        Output file path.
    email : str
        Entrez email address.
    db : str, default 'nucleotide'
        NCBI database name.
    start_acc : str, optional
        Accession from which fetching should resume.
    write_mode : str, default 'a'
        File mode used to open the output file.

    Returns
    -------
    None
    """
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

def filter_genbank_by_accession(input_file: str, output_file: str, accession_list: Iterable[str], verbose: bool = False) -> list[str]:
    """
    Filter a multi-entry GenBank file by accession.

    Parameters
    ----------
    input_file : str
        Path to the input GenBank file.
    output_file : str
        Path to the filtered output GenBank file.
    accession_list : list
        Versioned accession numbers to keep.
    verbose : bool, default False
        Whether progress messages should be printed.

    Returns
    -------
    list
        Requested accessions that were not found in the input file.
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

def append_genbank_from_gb(main_file: str, other_file: str, verbose: bool = False) -> None:
    """
    Append missing GenBank records from one file into another.

    Parameters
    ----------
    main_file : str
        Path to the GenBank file to update.
    other_file : str
        Path to the GenBank file providing candidate records.
    verbose : bool, default False
        Whether progress messages should be printed.

    Returns
    -------
    None
    """
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

def check_missing(gb_file: str,accession_list: Iterable[str]) -> Optional[list[str]]:
    """
    Check which accessions are absent from a GenBank file.

    Parameters
    ----------
    gb_file : str
        Path to the GenBank file to inspect.
    accession_list : iterable
        Accessions expected to be present.

    Returns
    -------
    list or None
        Missing accessions, or None when all requested records are present.
    """
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
def gb_to_fasta(input_file: str,output_file: str) -> None:
    """
    Convert a GenBank file to FASTA format.

    Parameters
    ----------
    input_file : str
        Input GenBank file path.
    output_file : str
        Output FASTA file path.

    Returns
    -------
    None
    """
    with open(output_file, "w") as fasta_file:
        SeqIO.write(SeqIO.parse(input_file, "genbank"), fasta_file, "fasta")    

def gb_feature_to_fasta(input_file: str, output_file: str, feature_type: str = 'CDS') -> None:
    """
    Export selected GenBank feature sequences to FASTA.

    Parameters
    ----------
    input_file : str
        Input GenBank file path.
    output_file : str
        Output FASTA file path.
    feature_type : str, default 'CDS'
        Feature type to extract.

    Returns
    -------
    None
    """
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

def get_gb(gb_file: str,accessions: Iterable[str]) -> list[Any]:
    """
    Retrieve selected GenBank records from a local file.

    Parameters
    ----------
    gb_file : str
        GenBank file path.
    accessions : iterable
        Record identifiers to extract.

    Returns
    -------
    list
        Matching sequence records.
    """
    outrecords=[]
    for record in SeqIO.parse(gb_file,'genbank'):
        if record.id in accessions:
            outrecords.append(record)
    return outrecords

def create_lookup(gbfile: str) -> set[str]:
    """
    Build a set of record identifiers from a GenBank file.

    Parameters
    ----------
    gbfile : str
        Path to the GenBank file.

    Returns
    -------
    set
        Set of GenBank record identifiers.
    """
    ids=set()
    for record in SeqIO.parse(gbfile,'genbank'):
        ids.add(record.id)
    return ids

def is_record(gbfile: str,accession: str) -> bool:
    """
    Check whether a GenBank file contains a record identifier.

    Parameters
    ----------
    gbfile : str
        GenBank file path.
    accession : str
        Record identifier to search for.

    Returns
    -------
    bool
        True when the record exists, otherwise False.
    """
    for record in SeqIO.parse(gbfile,'genbank'):
        if record.id == accession:
            return True
    return False
def fetch_genbank_list(access: list[str], email: str, db: str = 'nucleotide') -> list[Any]:
    """
    Fetch GenBank records from NCBI into memory.

    Parameters
    ----------
    access : list
        Accession numbers to fetch.
    email : str
        Entrez email address.
    db : str, default 'nucleotide'
        NCBI database name.

    Returns
    -------
    list
        Parsed GenBank sequence records.
    """
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

def write_gb(records: Iterable[Any],output: str) -> None:
    """
    Write GenBank records to a file.

    Parameters
    ----------
    records : iterable
        Sequence records to write.
    output : str
        Output GenBank file path.

    Returns
    -------
    None
    """
    with open(output,'w') as gbfile:
        SeqIO.write(records,gbfile,'genbank')

def append_genbank_from_list(main_file: str, records: Iterable[Any], verbose: bool = False) -> None:
    """
    Append missing GenBank records from an in-memory record list.

    Parameters
    ----------
    main_file : str
        Path to the GenBank file to update.
    records : iterable
        Sequence records to add when their identifiers are absent.
    verbose : bool, default False
        Whether progress messages should be printed.

    Returns
    -------
    None
    """
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

def fetch_fasta(access: list[str], file_name: str, email: str, db: str = 'nucleotide',start_acc: Optional[str] = None,write_mode: str = 'a') -> Optional[list[Any]]:
    """
    Fetch FASTA records from NCBI and append them to a file.

    Parameters
    ----------
    access : list
        Accession numbers to fetch.
    file_name : str
        Output file path.
    email : str
        Entrez email address.
    db : str, default 'nucleotide'
        NCBI database name.
    start_acc : str, optional
        Accession from which fetching should resume.
    write_mode : str, default 'a'
        File mode used to open the output file.

    Returns
    -------
    list or None
        Failed accession batches when an error occurs, otherwise None.
    """
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
def filter_fasta_by_accession(infile: str,access: Iterable[str],outfile: str) -> None:
    """
    Filter a FASTA file to records with selected accessions.

    Parameters
    ----------
    infile : str
        Input FASTA file path.
    access : iterable
        Accessions to retain.
    outfile : str
        Output FASTA file path.

    Returns
    -------
    None
    """
    outputs=[]
    with open(infile,'r') as handle:
        records=SeqIO.parse(handle,'fasta')
        for record in records:
            if record.id in access:
                outputs.append(record)
    with open(outfile,'w') as filtered:
        for elem in outputs:
            SeqIO.write(elem,filtered,'fasta')

def remove_records(fasta_file: str,accessions: Iterable[str],output_file: str) -> None:
    """
    Remove selected records from a FASTA file.

    Parameters
    ----------
    fasta_file : str
        Input FASTA file path.
    accessions : iterable
        Record identifiers to remove.
    output_file : str
        Output FASTA file path.

    Returns
    -------
    None
    """
    records=[]
    with open(fasta_file,'r') as handle:
        for record in SeqIO.parse(handle,'fasta'):
            if record.id not in accessions:
                records.append(record)
    with open(output_file,'w') as output:
        for record in records:
            SeqIO.write(record,output,'fasta')
