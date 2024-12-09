#!/bin/bash
help(){
    echo "#############################################################"
    echo "Install file for "
    echo "Requires ncbi-blast+, cd-hit, python 3, biopython and pandas"
    echo "------------------------------------------------------------"
    echo "To run: activate the appropriate conda fluid ENV"
    echo "./id_pipeline [DB_FASTA.fasta] [METADATA.csv]"
    echo "Creates dirtree and blast/cd-hit databases"
    echo "#############################################################"
}

FASTA=$1
METADATA=$2
ENV_NAME=fluid

if [ $# -ne 2 ]
then
    echo "ERROR: Invalid number of elements" >&2
    help
    exit 1
fi
if [[ "$CONDA_DEFAULT_ENV" != "$ENV_NAME"  ]]
then
    echo "ERROR: $ENV_NAME not active" >&2
	help
	exit 1
fi

echo "Creating project directories"

if [ ! -d "runs" ]
then
mkdir runs
fi
if [ ! -d "reports" ]
then
mkdir reports
fi
if [ ! -d "samples" ]
then
mkdir samples
fi
if [ ! -d "blast_db" ]
then
mkdir blast_db
fi
if [ ! -d "cluster_db" ]
then
mkdir cluster_db
fi
if [ ! -d "metadata" ]
then
mkdir metadata
fi

mv $METADATA metadata

echo "Creating BLAST database"

makeblastdb -in $FASTA -dbtype nucl -out ./blast_db/sequencesDnaInf

echo "Creating cd-hit cluster database"

cd-hit-est -i $FASTA -o ./cluster_db/infDNAClusters -c 0.99 -M 5000

echo "Creating cluster annotation file"

python3 clusterCharact.py

echo "Creating representative BLAST db"

makeblastdb -in ./cluster_db/infDNAClusters -dbtype nucl -out ./blast_db/infDNAClusters
