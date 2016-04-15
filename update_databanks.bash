#!/bin/bash

mkdir -p /data/fasta/ /data/blast/

build_templates () {

    python makeXrayPDBFinder2Fasta.py /data/fasta/templates.fa
    makeblastdb -in /data/fasta/templates.fa -dbtype prot -out /data/blast/templates
}

build_sprot () {

    rsync rsync.ebi.ac.uk::pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz /data/fasta/uniprot_sprot.fasta.gz
    gunzip -f /data/fasta/uniprot_sprot.fasta.gz

    makeblastdb -in /data/fasta/uniprot_sprot.fasta -dbtype prot -out /data/blast/uniprot_sprot
}

build_trembl () {

    rsync rsync.ebi.ac.uk::pub/databases/uniprot/knowledgebase/uniprot_trembl.fasta.gz /data/fasta/uniprot_trembl.fasta.gz
    gunzip -f /data/fasta/uniprot_trembl.fasta.gz

    makeblastdb -in /data/fasta/uniprot_trembl.fasta -dbtype prot -out /data/blast/uniprot_trembl
}

build_templates &
build_trembl &
build_sprot &

wait

/bin/echo -e "TITLE uniprot\nDBLIST /data/blast/uniprot_sprot /data/blast/uniprot_trembl" > /data/blast/uniprot.pal
