# Deployment

## Pre-requisites

* sudo apt-get install clustalw ncbi-blast+ rabbitmq-server python-pip

* install interproscan
  (see: https://code.google.com/p/interproscan/wiki/HowToDownload)
  and set its path in the config file under INTERPROSCAN

* install yasara
  (see: http://www.yasara.org/update.htm)
  and set the path of the yasara directory in the config file to YASARADIR

## Installation

* create the blast database of templates by first by running:
   makeXrayPDBFinder2Fasta.py and creating a fasta file, then run:
   makeblastdb -in <fasta> -dbtype prot -parse_seqids -out <database>
   finally, set the path of the database in the config file to TEMPLATESDB

* create the blast database of uniprot species by first running:
   uniprotSpecies.py and creating a fasta file, then run on each fasta:
   makeblastdb -in <fasta> -dbtype prot -parse_seqids -out <database dir>
   to create databases for each species
   finally, set the path of the database directory in the config file to
   SPECIESDBDIR

* copy all cmbi dssp files from ftp://ftp.cmbi.ru.nl/pub/molbio/data/dssp/
  and set the path to the dssp directory to DSSPDIR in the config file

* compile joanna lange's secondary structure-dependent alignment tool
  and set its path to MSA in the config file

* set MODELDIR, INTERPRODIR and EXECUTIONDIR to the directories to store
  models, interpro files and temporary runtime files respectively
