from flu_utils import seq_get,dict_to_fasta,pkl_save
import sys
def seq_enum(counter: int) -> str:
    """
    Generate a sequence identifier string in the format 'Seq_{counter}'.

    Parameters:
    counter (int): The counter value to be used in the sequence identifier.

    Returns:
    str: A sequence identifier string in the format 'Seq_{counter}'.
    """
    return f'>Seq_{counter}'

#Parsing arguments
if __name__ == '__main__':
    filename=sys.argv[1]
fasta_dict=seq_get(filename)

#creating a dictionary for header mappings
header_mapping={}
counter=1
for header in fasta_dict:
    key=seq_enum(counter)
    header_mapping[key]=header
    counter+=1
#converting to processed dictionary
out_dict={}
for key in header_mapping:
    out_dict[key]=fasta_dict[header_mapping[key]]

samples_path='samples/'
runs_path='runs/'
outname=filename.split('/')[-1]
outname=outname.replace('.fasta','')
outname='format_'+outname
dict_to_fasta(out_dict, f'{samples_path}{outname}')
pkl_save(header_mapping, f'{samples_path}{outname}_mappings')