#!/bin/bash
help(){
    echo "#############################################################"
    echo "Pipeline for identifying influenza genotypes and segments"
    echo "Requires ncbi-blast+, cd-hit, python 3, biopython and pandas"
    echo "------------------------------------------------------------"
    echo "To run: activate the appropriate conda ENV, containing cd-hit"
    echo "./id_pipeline [SAMPLE_NAME.fasta]"
    echo "Outputs reports to the reports/ directory"
    echo "#############################################################"
}

SECONDS=0 #timer

#PATHS
BLAST_DB=./blast_db
CLUSTER_DB=./cluster_db
METADATA=./metadata
REPORTS=./reports
SAMPLES_DIR=./samples
LOGS_DIR=./logs
RUN_DIR=./runs

#ARGS
SAMPLE=$1
OUTNAME=$(echo $SAMPLE | cut -d \. -f 1)
RUN_DT=$(date +"%Y-%m-%d_%H-%M-%S")
ID=0.99
CLUSTER_FASTA=infDNAClusters
CLUSTER_REFS=cluster_desc.txt
ENV_NAME=fluid
RM_PREV_RUN=true


#CREATING LOG AND ERROR LOG FILES
if [ ! -d $LOGS_DIR ]
then
    mkdir logs
fi
touch $LOGS_DIR/"$SAMPLE-$RUN_DT.stdout"
touch $LOGS_DIR/"$SAMPLE-$RUN_DT.stderr"
if [ ! -d $RUN_DIR ]
then
    mkdir runs
fi
#CLEARING PREVIOUS RUNS
if $RM_PREV_RUN;
then
	find $RUN_DIR -type f ! -name "*.pkl" -exec rm -f {} +
    find $RUN_DIR -type f -name "*.pkl" -exec rm -f {} +
    find $SAMPLE_DIR -type f -name "to_blast.fasta" -exec rm -f {} +
    find $SAMPLE_DIR -type f -name "to_reblast.fasta" -exec rm -f {} +
    find $REPORTS -type f -name "to_blast_report.txt" -exec rm -f {} +
    find $REPORTS -type f -name "to_reblast_report.txt" -exec rm -f {} +
fi

#PATH AND FILE CHECKS
if [ ! -d $METADATA ]
then
    echo "ERROR: $METADATA dir not found" >&2
    help
    exit 1
fi
if [ ! -d $SAMPLES_DIR ]
then
    echo "ERROR: $SAMPLES_DIR dir not found" >&2
    help
    exit 1
fi
if [ ! -d $BLAST_DB ]
then
    echo "ERROR: $BLAST_DB dir not found" >&2
    help
    exit 1
fi
if [ ! -d $CLUSTER_DB ]
then
    echo "ERROR: $CLUSTER_DB dir not found" >&2
    help
    exit 1
fi
if [ ! -d $REPORTS ]
then
    echo "ERROR: $REPORTS dir not found" >&2
    help
    exit 1
fi
if [ ! -e $SAMPLES_DIR/$SAMPLE ]
then
    echo "ERROR: $SAMPLE file not found" >&2
    help
    exit 1
fi
if [ $# -ne 1 ]
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
################################
#CD-HIT RUN
cd-hit-est-2d -i $CLUSTER_DB/$CLUSTER_FASTA -i2 $SAMPLES_DIR/$SAMPLE -o $RUN_DIR/$OUTNAME -c $ID >>$LOGS_DIR/"$SAMPLE-$RUN_DT.stdout"\
 2>>$LOGS_DIR/"$SAMPLE-$RUN_DT.stderr"

#RETURNS ASSIGNMENTS FROM CD-HIT
python3 clusterAssign.py $RUN_DIR/"$OUTNAME.clstr" $SAMPLES_DIR/$SAMPLE >>$LOGS_DIR/"$SAMPLE-$RUN_DT.stdout"\
 2>>$LOGS_DIR/"$SAMPLE-$RUN_DT.stderr"


#GENERATES $SAMPLE_clust.txt report with: sample; ID; Assign_Cluster; Cluster_rep; genotype; segment; host
python3 clusterCompile.py $RUN_DIR/dict_sample.pkl $RUN_DIR/"$OUTNAME.assign" >>$LOGS_DIR/"$SAMPLE-$RUN_DT.stdout"\
 2>>$LOGS_DIR/"$SAMPLE-$RUN_DT.stderr"

#MINES THE $SAMPLE_clust.txt REPORT FOR UNASSIGNED SAMPLES

python3 clusterMiner.py "$OUTNAME""_clust.txt" $SAMPLE >>$LOGS_DIR/"$SAMPLE-$RUN_DT.stdout"\
 2>>$LOGS_DIR/"$SAMPLE-$RUN_DT.stderr"

#SINGLE-MATCH BLAST SCRIPT ON UNASSIGNED SAMPLES
if [ -s "$SAMPLES_DIR/to_blast.fasta" ]
then
    python3 blastCluster.py to_blast.fasta >>$LOGS_DIR/"$SAMPLE-$RUN_DT.stdout" 2>>$LOGS_DIR/"$SAMPLE-$RUN_DT.stderr"
    if [ -s "$SAMPLES_DIR/to_reblast.fasta" ]
    then 
        python3 bestBlast2.py to_reblast.fasta  >>$LOGS_DIR/"$SAMPLE-$RUN_DT.stdout" 2>>$LOGS_DIR/"$SAMPLE-$RUN_DT.stderr"
    fi

fi

#CONFORMING DATA IN FINAL REPORT
python3 reportGenerator.py to_blast_bclust_report.txt to_reblast_report.txt $OUTNAME 2>>$LOGS_DIR/"$SAMPLE-$RUN_DT.stderr"

RUNTIME=$SECONDS
echo "Script ended in $RUNTIME seconds"
