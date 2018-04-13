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
MODELS_BLAST=$BLAST_DIR/models

build_models () {

    $PYTHON make_models_fasta.py $MODELS_FASTA
    $MAKEBLASTDB -in $MODELS_FASTA -dbtype prot -out $MODELS_BLAST
}

TEMPLATES_FASTA=$FASTA_DIR/templates.fa
TEMPLATES_BLAST=$BLAST_DIR/templates

build_templates () {

    $PYTHON make_templates_fasta.py $TEMPLATES_FASTA
    $MAKEBLASTDB -in $TEMPLATES_FASTA -dbtype prot -out $TEMPLATES_BLAST
}

SPROT_FASTA=$FASTA_DIR/uniprot_sprot.fasta
SPROT_BLAST=$BLAST_DIR/uniprot_sprot

build_sprot () {

    $RSYNC rsync.ebi.ac.uk::pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz $SPROT_FASTA.gz
    gunzip -f $SPROT_FASTA.gz

    # To prevent warnings, remove all titles from the fasta.
    sed -i 's/^>\([^ ]\+\) .*$/>\1/' $SPROT_FASTA

    $MAKEBLASTDB -in $SPROT_FASTA -dbtype prot -out $SPROT_BLAST
}

TREMBL_FASTA=$FASTA_DIR/uniprot_trembl.fasta
TREMBL_BLAST=$BLAST_DIR/uniprot_trembl

build_trembl () {

    $RSYNC rsync.ebi.ac.uk::pub/databases/uniprot/knowledgebase/uniprot_trembl.fasta.gz $TREMBL_FASTA.gz
    gunzip -f $TREMBL_FASTA.gz

    # To prevent warnings, remove all titles from the fasta.
    sed -i 's/^>\([^ ]\+\) .*$/>\1/' $TREMBL_FASTA

    $MAKEBLASTDB -in $TREMBL_FASTA -dbtype prot -out $TREMBL_BLAST
}

build_models &
build_templates &
build_trembl &
build_sprot &

wait

/bin/echo -e "TITLE uniprot\nDBLIST $SPROT_BLAST $TREMBL_BLAST" > $BLAST_DIR/uniprot.pal
