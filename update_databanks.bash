#!/bin/bash

PYTHON=/usr/local/bin/python
MAKEBLASTDB=/usr/bin/makeblastdb
RSYNC=/usr/bin/rsync

DATA_DIR=/data
FASTA_DIR=$DATA_DIR/fasta
BLAST_DIR=$DATA_DIR/blast
MODEL_DIR=$DATA_DIR/models

mkdir -p $FASTA_DIR $BLAST_DIR

MODELS_FASTA=$FASTA_DIR/models.fa

build_models () {

    mkdir -p /tmp/models-blast

    $PYTHON make_models_fasta.py $MODELS_FASTA
    $MAKEBLASTDB -in $MODELS_FASTA -dbtype prot -out /tmp/models-blast/models
    mv -f /tmp/models-blast/models.* $BLAST_DIR/
}

TEMPLATES_FASTA=$FASTA_DIR/templates.fa

build_templates () {

    mkdir -p /tmp/templates-blast

    $PYTHON make_templates_fasta.py $TEMPLATES_FASTA
    $MAKEBLASTDB -in $TEMPLATES_FASTA -dbtype prot -out /tmp/templates-blast/templates
    mv -f /tmp/templates-blast/templates.* $BLAST_DIR/
}

SPROT_FASTA=$FASTA_DIR/uniprot_sprot.fasta

build_sprot () {

    $RSYNC rsync.ebi.ac.uk::pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz $SPROT_FASTA.gz
    gunzip -f $SPROT_FASTA.gz

    # To prevent warnings, remove all titles from the fasta.
    sed -i 's/^>\([^ ]\+\) .*$/>\1/' $SPROT_FASTA

    mkdir -p /tmp/sprot-blast

    $MAKEBLASTDB -in $SPROT_FASTA -dbtype prot -out /tmp/sprot-blast/uniprot_sprot
    mv -f /tmp/sprot-blast/uniprot_sprot.* $BLAST_DIR/
}

TREMBL_FASTA=$FASTA_DIR/uniprot_trembl.fasta

build_trembl () {

    $RSYNC rsync.ebi.ac.uk::pub/databases/uniprot/knowledgebase/uniprot_trembl.fasta.gz $TREMBL_FASTA.gz
    gunzip -f $TREMBL_FASTA.gz

    # To prevent warnings, remove all titles from the fasta.
    sed -i 's/^>\([^ ]\+\) .*$/>\1/' $TREMBL_FASTA

    mkdir -p /tmp/trembl-blast

    $MAKEBLASTDB -in $TREMBL_FASTA -dbtype prot -out /tmp/trembl-blast/uniprot_trembl
    mv -f /tmp/trembl-blast/uniprot_trembl.* $BLAST_DIR/
}

build_models &
build_templates &
build_trembl &
build_sprot &

wait

/bin/echo -e "TITLE uniprot\nDBLIST $BLAST_DIR/uniprot_sprot $BLAST_DIR/uniprot_trembl" > $BLAST_DIR/uniprot.pal
