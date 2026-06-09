# AFluID
Pipeline for Influenza A virus segment and genotype identification though homolgy clustering and local alignment

## Installation

### Path policy
- Run maintenance commands from the AFluID project root. This includes `install.py` and `update.py`.
- The configuration file may be passed as either a relative or absolute path, but AFluID resolves it to an absolute path at startup and uses that absolute path for the rest of the run.
- Paths defined under `[Paths]` in `config.ini` may be absolute or relative.
- Absolute paths under `[Paths]` are used as written.
- Relative paths under `[Paths]` are resolved relative to the project root.
- Maintenance input files required by `install.py` and `update.py` remain rooted at the project root.

### Requirements
1. Anaconda or Miniconda installed and wsl or unix/linux command line
2. A fasta file with several fully assembled sequences of multiple influenza segments: I got the dataset I used from NCBI virus database.
3. A metadata .csv file describing each acession in **at a minimum** segment, genotype and host. I have a manually curated version in this repository (ver. 4) that has the header fields necessary for the python indexing to work. The dataset used for this work is available [here](https://zenodo.org/records/20431631).

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

   **Warning:** the AFluID database is built and updated from complete assemblies with correctly annotated sequences only. Assemblies that are incomplete, contain low-quality sequence, or have missing or incorrect metadata should be rejected before building or updating the database, otherwise cluster annotations and downstream reports may become unreliable.

4. Use **dbSetup.py** to filter out all the sequences that are not included within the metadata table:
```
python3 dbSetup.py  [PATH_TO_METADATA_FILE] [PATH_TO_FASTA_FILE] [OUTPUT_FILE_NAME]
``` 
Place the maintenance input files in the project root and keep their names aligned with the entries in **config.ini**.
  * Should you have all the sequences you need the script will output a file with all the fasta sequences that you have that are on the metadata file
  * If your sequence FASTA file is incomplete, the script will output a .txt file with the missing sequences' accession numbers for you to download and concatenate to your main file

5. If you haven't done so yet, open the config.ini file in a notepad program (notepad, notepad++, sublime text...) or nano to change any settings you deem necessary for your analysis. **You also need to specify your email and conda environment path**.

   Required maintenance input files in the project root:
   * metadata CSV defined by `Filenames.metadata`
   * reference GenBank file defined by `Filenames.ref_db`
   * `geo_tree.json`
   * `host_taxonomy_tree.json`

6. Run the installation script using the command:
```
python3 install.py config.ini
```
   `install.py` must be run from the project root. Paths defined under `[Paths]` in `config.ini` may be absolute or relative: absolute paths are used as written, while relative paths are resolved relative to the project root.
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
* Single-sample (*argument: -ss/--single_sample (on/off) - default:on*):
    Single-sample is a sub-mode of consensus mode used for when a user has a fasta file with 1-8 sequences of one sample. This mode will allow for more quality checks on the final report. For multi-sample fasta files please turn this option off.
* Turn off tools (*argument: -off/--turn_off (nextclade,getref,flumut)*):
    This argument allows the user to turn off specific tools thar could be turned on by default.
* Remove previous (*argument: -rm/--remove_previous (on/off) - default:on*):
    The pipline profuces intermediate outputs on the ```runs/``` and ```samples/``` paths. this will clear those files in the beginning of each run. If you are running batch or parallell analysis, turn off this feature.
* Output directory (*argument: --outdir PATH*):
     Report outputs and exported reference files can be redirected to a custom root directory. Relative paths are resolved from the current working directory; absolute paths are accepted as-is.
* Output name (*argument: -o/--output_name NAME*):
     In direct single-file runs, this optionally overrides the default output stem used for the output folder and generated report filenames in both ```-ss on``` and ```-ss off``` modes. The name is used exactly as provided, without timestamps or hex suffixes, so an existing target raises an error unless ```--replace``` is used. This option cannot be used with ```--batch``` or ```--batch_fasta```.
* Batch name (*argument: -bn/--batch_name NAME*):
     In batch mode or successful batch-fasta demultiplexing, this defines the batch output directory name. If omitted, AFluID uses the current timestamped batch naming behaviour. If ```--batch_fasta``` falls back to multi-sample mode through ```--force-bf```, ```--batch_name``` is ignored.
* Batch FASTA demultiplexing (*arguments: -bf/--batch_fasta, --bf-regex REGEX, --force-bf*):
     This optional mode allows AFluID to inspect one multifasta input, gather records by a repeated sample identifier captured by a user-provided regex, and run the result as batch single-sample mode. ```--bf-regex``` is mandatory with ```--batch_fasta``` and must capture only the sample identifier, either with a named group ``sample`` or with the first capture group. If every header matches, the file is demultiplexed and processed as ```-b -ss on```. If any header fails, the run stops unless ```--force-bf``` is provided, in which case the original file is processed as a multi-sample FASTA (effective ```-ss off```). Original FASTA headers are preserved unchanged in the demultiplexed files. Demultiplexed FASTA filenames are sanitized to filesystem-safe names derived from the captured sample identifier, and a short deterministic suffix is added only if two captured identifiers would otherwise collide on disk.
* Replace outputs (*argument: --replace*):
     Prevent accidental overwrites by default. If the target output directory already exists, the run will fail unless ```--replace``` is provided, in which case the existing output directory is deleted and recreated before the run starts.

The main pipeline may be launched from outside the project root, but any relative paths defined under `[Paths]` in `config.ini` are still interpreted relative to the project root, not the shell working directory.

To run the pipeline **place the sample fasta files in the ```samples/``` folder**. You only need to specify the name of the file to be analysed as the pipeline will automatically search for the sample files in the ```samples/``` folder

You can find intermediate output files in run-specific temporary workspaces under ```runs/``` while the pipeline is executing, and logs in the ```logs/``` directory. Successful runs clean up their temporary ```runs/``` workspace automatically.

In single-sample mode, reports are written to per-sample folders under the configured ```reports/``` directory and exported references are written to per-sample folders under the configured ```references/``` directory. When ```--outdir``` is used, reports are written to ```<outdir>/<sample>/reports/``` and exported references to ```<outdir>/<sample>/references/```.

In batch mode, reports and exported references are grouped under a batch directory. With the default paths, reports are written under ```<reports>/<batch_name>/<sample>/``` and exported references under ```<references>/<batch_name>/<sample>/```. When ```--outdir``` is used, reports are written to ```<outdir>/<batch_name>/<sample>/reports/``` and exported references to ```<outdir>/<batch_name>/<sample>/references/```. The batch summary, batch index, and batch ZIP are written at the batch root, and the HTML index and ZIP preserve the nested sample folder structure.

Sample usage (one sample per file, default config file, single sample mode, no forces, no offs, consensus mode):

```
python3 main.py -f  SAMPLE_FASTA.fasta -m consensus
```

Sample usage (direct run with a custom output stem and collision checking):

```bash
python3 main.py -f SAMPLE_FASTA.fasta -m consensus -ss on -o SAMPLE_ALIAS
```

Sample usage (demultiplex a multifasta into a temporary batch using a mandatory regex that captures only the sample identifier):

```bash
python3 main.py -f MY_MULTIFASTA.fasta -m consensus -bf --bf-regex '(?P<sample>EPI_ISL_\d+)'
```

Example using the first capture group instead of a named group:

```bash
python3 main.py -f MY_MULTIFASTA.fasta -m consensus -bf --bf-regex '(A/[^_ ]+/[A-Za-z0-9-]+/\d{4})'
```

Example forcing a non-demultipliable multifasta to run as multi-sample mode:

```bash
python3 main.py -f MY_MULTIFASTA.fasta -m consensus -bf --bf-regex '(?P<sample>[^|]+)' --force-bf
```

## Pipeline steps
This pipeline comprises several steps:
1. The pipeline preprocesses the fasta file in order to remove large headers and metacharacters. Will create a mapping dictionary to map the conformed headers against the original ones. This step will also filter sequences below and above preconfigured length thresholds. Every discarded sequence will be stored in the flags dictionary.

2. The sequences will be subjected to CD-HIT clustering

3. The CD-HIT outputs are parsed and a cluster report is outputted.

4. The unassigned sequences will be subjected to a BLAST step against cluster representatives.

5. Every unassigned samples for this step will be subjected to a BLAST step against the bulk local database.

6. The stored data form the previous reports will then be compiled into the final ID report.

7. The pipeline will then parse the report and redirect every sequence to the posterior analysis (flumut, genin, nextclade or get_references)

The following diagram illustrates the pipeline workflow:
[Pipeline workflow (light)](Diagrama_pipeline.drawio.png) or [Pipeline workflow (dark)](Diagrama_pipeline.drawio.dark.png)

## Mutations of interest
The pipeline parses and reports the following mutations of interest. This reference table is derived from the current hard-coded mutation set in `structures.py`, specifically `muts_interest` and `muts_loci_meaning`.

| Segment/Locus | Mutation | Biological effect |
| --- | --- | --- |
| PB2 | 9N | Observed in domestic poultry, Uganda; increased virulence in mice |
| PB2 | 199S | Increased virulence in mice |
| PB2 | 271A | Increase polymerase activity in the presence of avian acidic nuclear phosphoprotein 32 (ANP32) proteins; increase mortality and airborne transmission in ferrets; increased polymerase activity in avian and mammalian cell line |
| PB2 | 292V | Observed in H3N8 samples from humans and birds; increased polymerase activity in a mammalian cell line, increased virulence in mice |
| PB2 | 526R | Increased polymerase activity in mammalian cell line and virulence in mice |
| PB2 | 588I | Enhances 2009 H1N1 pandemic influenza virus virulence by increasing viral replication and exacerbating PB2 inhibition of beta-interferon expression |
| PB2 | 588V | Increased polymerase activity and replication in mammalian and avian cell lines, increased virulence in mice |
| PB2 | 591K | Increased polymerase activity in mammalian and avian cell line, increased replication in a mammalian cell line, increased virulence in mice; observed in human isolates |
| PB2 | 591R | Enhanced mammalian ANP32 adaptation (non-biased) |
| PB2 | 627K | Observed in human isolate; increased polymerase activity and replication in mammalian cell line, increased virulence in mice and ferrets; contributes to airborne transmission of influenza A viruses (IAVs) in ferrets and contact transmission in guinea pigs; mammalian ANP32-specific adaptation (ANP32B biased) |
| PB2 | 627V | Observed in human isolate; mammalian ANP32-specific adaptation (ANP32B biased); observed in human samples in Shenzhen; increased polymerase activity and replication in mammalian cell lines, increased virulence in mice |
| PB2 | 631L | Most dominant mutation in mouse-adapted virus that strongly upregulated viral polymerase activity and played a critical role in the enhancement of virus replication and disease severity in mice; observed in dairy cattle USA 2024 |
| PB2 | 701N | Increased polymerase activity, enhanced replication efficiency, increased virulence and contact transmission in guinea pigs; increased polymerase activity in mammalian cell line; increased viral replication and virulence in mice; mammalian ANP32-specific adaptation (non-biased); observed in seals |
| PB2 | 702N | Observed in human isolate; increased polymerase activity in human and avian cells |
| PB2 | 740N | Increase polymerase activity in the presence of avian ANP32 proteins, in human and avian cells |
| PB1-F2 | 66S | Increased virulence, replication efficiency and antiviral response in mice |
| PA | 38T | Resistance to antivirals: Baloxavir Marboxil |
| PA | 85I | Enables guanine-rich sequence binding factor (GRSF1) to enhance the cytosolic accumulation and translation of a subset of viral mRNAs; enhanced replication in human A549 cells |
| PA | 97I | Increased polymerase activity in mammalian cell line and enhanced replication and virulence in mice; increased polymerase activity in mammalian |
| PA | 186S | Enables GRSF1 to enhance the cytosolic accumulation and translation of a subset of viral mRNAs |
| PA | 336M | Significantly enhanced pathogenicity in a mouse model |
| PA | 356R | Human host marker; Increase polymerase activity and enhanced replication in a mammalian cell line, increased virulence in mice |
| PA | 552S | Enhanced viral RNA-dependent RNA polymerase (vRdRp) activity and viral replication in vitro |
| HA | 156A | Increased virus binding to α2,6, increased transmission in guinea pigs |
| HA | 156V | Observed in seal samples; decreased virulence in mice and increased affinity for the human-type receptor |
| HA | 186D | In combination with HA:221D; Switch in the receptor specificity from avian-type to human-type receptor |
| HA | 186V | Enhances binding affinity to mammalian cells and replication in mammalian cells; enhanced replication in mice |
| HA | 208T | Increased the viral replication in avian and mammalian cells; enhance viral replication in mice |
| HA | 221D | Increased virus binding to α2,6 |
| HA | 222L | Increased virus binding to α2,6; transmitted via aerosol among guinea pigs; enhanced replication in mammalian cells and ferrets, enhanced contact transmission in ferrets; Loss of binding to α2,3; Increased acid and thermal stability |
| HA | 224S | Increased binding to α2,6; increased viral replication in mammalian cells and virulence in mice; observed in human isolate; increase in the ability of the virus to infect mammals; decreased virus binding to α2,3 |
| NP | 52H | Role in Butyrophilin Subfamily 3 Member A3 (BTN3A3) evasion |
| NP | 52N | Significant mutation observed also in H5N1; role in BTN3A3 evasion |
| NP | 100I | Resistance against Myxovirus resistance protein A (MxA) |
| NP | 100V | Resistance against Myxovirus resistance protein A (MxA) |
| NP | 283P | Resistance against Myxovirus resistance protein A (MxA) |
| NP | 313V | Role in Butyrophilin Subfamily 3 Member A3 (BTN3A3) evasion |
| NP | 313Y | Role in Butyrophilin Subfamily 3 Member A3 (BTN3A3) evasion |
| NA | 119V | Resistance to antivirals: Zanamivir |
| NA | 223V | Resistance to antivirals: Laninamivir |
| NA | 275Y | Resistance to antivirals: Oseltamivir and Peramivir |
| NA | 292K | Resistance to antivirals: Zanamivir |
| NA | 399R | Mutation of 2SBS that are detrimental to the cleavage of sialosides linked to fetuin or transferrin; observed in seals |
| NA | 432E | Decreased cleavage of fetuin-containing α2,3-linked sialic acids (SIAs) but not that of monovalent substrates or of transferrin containing only α2,6-linked SIAs |
| M1 | 95K | Resistant to TRIM21, increases replication and pathogenicity in mice |
| M2 | 27A | Increased resistance to adamantanes |
| M2 | 31N | Increased resistance to adamantanes |

If `muts_interest` or `muts_loci_meaning` changes in `structures.py`, this table should be updated to stay aligned with the mutations currently parsed by AFluID.

## Update
The user can and should update AFluID to better represent the genomic diversity of the current flu season or the new circulating variants of non-seasonal Influenza A.
The update can be made using the update.py module. `update.py` is a maintenance command and must be run from the project root.
The user should go to NCBI Virus and download the sequence metadata with matching headers to the metadata, and keeping the column header case as is. The values should be from the date os the provious updato to the present. A fasta file with the same sequences should also be downloaded, from the same database.

**Warning:** update inputs should follow the same curation rules as the original database build. Only complete assemblies with correctly annotated sequences should be added. Incomplete assemblies, low-quality sequence, or records with missing or incorrect metadata must be excluded from the update set.

These files will need to be placed in the 
``update/`` directory and declared within the command run:
```
python3 update.py -c CONFIG_FILE -ff UPDATE_FASTA -mf UPDATE_CSV
```
As with installation, the configuration file is normalized to an absolute path at startup. Paths defined under `[Paths]` in `config.ini` may be absolute or relative: absolute paths are used as written, and relative paths are resolved relative to the project root.
Depending on the size of the sequence database the update can take a long while so it's recommended you use a screen if working on a shared server.

Update validation reports are written to the dedicated path configured as ```update_reports``` under ```[Paths]``` in ```config.ini```, so update runs stay isolated from normal pipeline report outputs. Older configs that do not define ```update_reports``` fall back to ```update/reports``` automatically.
