# AFluID
Pipeline for Influenza A virus segment and genotype identification though homolgy clustering and local alignment

## Installation

### Requirements
1. Anaconda or Miniconda installed and wsl or unix/linux command line
2. A fasta file with several fully assembled sequences of multiple influenza segments: I got the dataset I used from NCBI8
3. A metadata .csv file describing each acession in **at a minimum** segment, genotype and host. I have a manually curated version in this repository (ver. 4)

Begin by cloning this repository using the command:

```
git clone AFluID
```

Then set up a conda environment using the command
```
conda env create -f env.yml
```
