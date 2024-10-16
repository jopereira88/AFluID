# AFluID
Pipeline for Influenza A virus segment and genotype identification though homolgy clustering and local alignment

## Installation

### Requirements
1. Anaconda or Miniconda installed and wsl or unix/linux command line
2. A fasta file with several fully assembled sequences of multiple influenza segments: I got the dataset I used from NCBI virus database.
3. A metadata .csv file describing each acession in **at a minimum** segment, genotype and host. I have a manually curated version in this repository (ver. 4) that has the header fields necessary for the python indexing to work

### Procedure
1. Begin by cloning this repository using the command:
```
git clone https://github.com/jopereira88/AFluID.git
```
2. Then set up a conda environment using the command
```
conda env create -f env.yaml
```
3. Run the installation script install.sh using the following commands, don't forget to decompress the metadata file, if you use the one provided.
```
chmod +x install.sh     #do this the first time to make the file executable
./install.sh YOUR_FASTA_FILE.fasta METADATA_FILE.csv
``` 
4. This will create the directory tree and the database files. **Depending on the size of you data this step can take a fair bit of time**

## Running the pipeline

This pipeline will identify nucleotide sequences of influenza A into segments and genotypes using the two following methods:

1. The samples will be assigned into clusters using cd-hit-est-2d with an identity thresshold of 99%. These clusters are annotated for segment, genoypes and hosts, therefore metadata completeness is paramount to the reports. The samples that are unassigned to clusters will then be subjected to BLAST against a local database.
2. BLAST will output the best match, identity and the results will be annotated using sequence metadata. In this case, only one segment, genotype and host will be outputted if available.
3. Any unassigned samples by BLAST will be outputted so the user is able to perform online BLAST.

To run the pipeline place the sample fasta files in the ```samples/``` directory and run these commands:
```
chmod +x id_pipeline.sh      #do this the first time to make the file executable
./id_pipeline.sh SAMPLE_FASTA_NAME.fasta
```

You can find intermediate output files in the ```runs/``` directory, logs in the ```logs/``` directory and the intermediate and final reports in the ```reports/``` directory.

## Pipeline steps

