# AFluID
Pipeline for Influenza A virus segment and genotype identification though homolgy clustering and local alignment

## Installation

### Requirements
1. Anaconda or Miniconda installed and wsl or unix/linux command line
2. A fasta file with several fully assembled sequences of multiple influenza segments: I got the dataset I used from NCBI virus database.
3. A metadata .csv file describing each acession in **at a minimum** segment, genotype and host. I have a manually curated version in this repository (ver. 4) that has the header fields necessary for the python indexing to work
4. A genbank .gb file that can be empty (just use the command ```touch YOUR_GB_REF_FILE.gb```)

### Procedure
1. Begin by cloning this repository using the command:
```
git clone https://github.com/jopereira88/AFluID.git
```
2. Then set up a conda environment using the command
```
conda env create -f env.yaml
```
3. If you haven't downloaded the dataset yet I recommend using NCBI virus, searching for **Influenza A virus**, filter for **Complete Assemblies** and download all sequences. Also don't forget to **Activate the conda environment using:** ```conda activate fluid```

4. Use **dbSetup.py** to filter out all the sequences that are not included within the metadata table:
```
python3 dbSetup.py  [PATH_TO_METADATA_FILE] [PATH_TO_FASTA_FILE] [OUTPUT_FILE_NAME]
``` 
I recommend you use the cwd as a path for now and the names of metadata and output files equal to the ones in **config.ini**
  * Should you have all the sequences you need the script will output a file with all the fasta sequences that you have that are on the metadata file
  * If your sequence FASTA file is incomplete, the script will output a .txt file with the missing sequences' accession numbers for you to download and concatenate to your main file

5. If you haven't done so yet you should open the config.ini file in a notepad program (notepad, notepad++, sublime text...) or nano to change any settings you deem necessary for your analysis.

6. Run the installation script using the command:
```
python3 install.py config.ini
```
7. This will create the directory tree and the database files. **Depending on the size of your data and the speed of your computer every step can take a fair bit of time**

## Running the pipeline

This pipeline will identify nucleotide sequences of influenza A into segments and genotypes/subtypes using the two following methods:

1. The samples will be assigned into clusters using cd-hit-est-2d with an identity thresshold of 99%. These clusters are annotated for segment, genoypes and hosts, therefore metadata completeness is paramount to the reports. The samples that are unassigned to clusters will then be subjected to BLAST against a local database.
2. BLAST will output the best match, identity and the results will be annotated using sequence metadata. In this case, only one segment, genotype and host will be outputted if available.
3. Any unassigned samples by BLAST will be outputted so the user is able to perform online BLAST.
4. The pipeline works with an argparse argument calling system, there are 2 mandatory arguments: *-f/--filename (your_fasta_file_name.fasta)* and *-m/--mode*

The pipeline also runs 2 different modes and its behaviour can be altered by several different arguments:
* Contig mode (_argument: -m/--mode contig_):
    This mode should be used on _drafted contig files_ and will filter sequences by size, identify them and then get the genbank (.gb) and fasta (.fasta) files of the closest reference (using the getReference() function) in order to build assembly consensus.
* Consensus mode (_argument: -m/--mode consensus_):
    This mode will should be used for consensus sequences and will filter sequences by size, identify them and run several tools: flumut on any sequence identified as H5Nx; Nextclade on every H1,H3 and H5 HA sequences.
* Forcing tools (_argument: -ff/--force (flumut,getref)_):
    The argument force will allow the user to perform additional bioinformatics analysis that are not offered as defaults by the modes, for example, it will allow flumut analysis on drafted contigs or getref on consensus samples. Nextclade is not forceable at this moment.
* Single-sample (*argument: -ss/-single_sample (on/off) - default:on*):
    Single-sample is a sub-mode of consensus mode used for when a user has a fasta file with 1-8 sequences of one sample. This mode will allow for more quality checks on the final report. For multi-sample fasta files please turn this option off.
* Turn off tools (*argument: -off/--turn_off (nextclade,getref,flumut)*):
    This argument allows the user to turn off specific tools thar could be turned on by default.
* Remove previous (*argument: -rm/--remove_previous (on/off) - default:on*):
    The pipline profuces intermediate outputs on the ```runs/``` and ```samples/``` paths. this will clear those files in the beginning of each run. If you are running batch or parallell analysis, turn off this feature.

To run the pipeline place the sample fasta files in the ```samples/``` directory and run these commands:

You can find intermediate output files in the ```runs/``` directory, logs in the ```logs/``` directory and the intermediate and final reports in the ```reports/``` directory.

## Pipeline steps

