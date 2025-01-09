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
min_seq_len=200
max_seq_len=3500
verbose=True

#removing length outliers
to_pop=[]
for key in fasta_dict:
    if len(fasta_dict[key]) < min_seq_len or len(fasta_dict[key]) > 2380:
        to_pop.append(key)
for item in to_pop:
    fasta_dict.pop(item)
if verbose:
    if len(to_pop)==1:
        print(f'{len(to_pop)} sequence was removed for being outside the length thresholds:')
        #print(to_pop[0])
    elif to_pop==[]:
        print('All sequences were within length thresholds')
    else:
        print(f'{len(to_pop)} sequences were removed for being outside the length thresholds:')
        #print(to_pop)

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
outname_fasta='format_'+outname
dict_to_fasta(out_dict, f'{samples_path}{outname_fasta}')
pkl_save(header_mapping, f'{samples_path}{outname}_mappings')