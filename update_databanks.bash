#!/bin/bash

mkdir -p /data/fasta/ /data/blast/

# Only copy this file the first time, don't overwrite.
if ! [ -f /data/blacklisted_templates ] ; then
    cp blacklisted_templates /data/blacklisted_templates
fi

build_models () {

    python models-fasta.py /data/models /mnt/chelonium/hg/models /data/fasta/models.fa
    makeblastdb -in /data/fasta/models.fa -dbtype prot -out /data/blast/models
}

build_templates () {

    python makeXrayPDBFinder2Fasta.py /data/fasta/templates.fa
    makeblastdb -in /data/fasta/templates.fa -dbtype prot -out /data/blast/templates
}

build_sprot () {

    rsync rsync.ebi.ac.uk::pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz /data/fasta/uniprot_sprot.fasta.gz
    gunzip -f /data/fasta/uniprot_sprot.fasta.gz

    # To prevent warnings ..
    sed -i 's/^>\([^ ]\+\) .*$/>\1/' /data/fasta/uniprot_sprot.fasta

    makeblastdb -in /data/fasta/uniprot_sprot.fasta -dbtype prot -out /data/blast/uniprot_sprot
}

build_trembl () {

    rsync rsync.ebi.ac.uk::pub/databases/uniprot/knowledgebase/uniprot_trembl.fasta.gz /data/fasta/uniprot_trembl.fasta.gz
    gunzip -f /data/fasta/uniprot_trembl.fasta.gz

    # To prevent warnings ..
    sed -i 's/^>\([^ ]\+\) .*$/>\1/' /data/fasta/uniprot_trembl.fasta

    makeblastdb -in /data/fasta/uniprot_trembl.fasta -dbtype prot -out /data/blast/uniprot_trembl
}

build_models &
build_templates &
build_trembl &
build_sprot &

wait

/bin/echo -e "TITLE uniprot\nDBLIST /data/blast/uniprot_sprot /data/blast/uniprot_trembl" > /data/blast/uniprot.pal
