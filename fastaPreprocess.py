from flu_utils import seq_get,dict_to_fasta
import sys
if __name__ == '__main__':
    filename=sys.argv[1]
fasta_dict=seq_get(filename)
samples_path='samples/'
outname=filename.split('/')[-1]
outname=outname.replace('.fasta','')
outname='format_'+outname
dict_to_fasta(fasta_dict, f'{samples_path}{outname}')