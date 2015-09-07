#!/usr/bin/python

# This is a script for building models without using celery or the restful service.
# It can be run from console and uses the same modeling procedure as the restful service does.

import logging
import sys
import os
import traceback

sys.path.append (os.path.dirname (os.path.abspath (__file__)))
from hommod_rest.services.model import modeler
from hommod_rest.services.modelutils import parseFasta
import hommod_rest.default_settings as config

_log = logging.getLogger(__name__)


def init():

    # Initialise services:

    from hommod_rest.services.model import modeler
    modeler.yasara_dir = config.YASARADIR
    modeler.execution_root_dir = config.EXECUTIONDIR
    modeler.model_root_dir = config.MODELDIR

    from hommod_rest.services.secstr import secstr
    secstr.dssp_dir = config.DSSPDIR
    secstr.yasara_dir = config.YASARADIR

    from hommod_rest.services.interpro import interpro
    interpro.interproExe = config.INTERPROSCAN
    interpro.storageDir = config.INTERPRODIR

    from hommod_rest.services.blast import blaster
    blaster.templatesDB = config.TEMPLATESDB
    blaster.uniprotspeciesDir = config.SPECIESDBDIR
    blaster.blastpExe = config.BLASTP

    from hommod_rest.services.align import aligner
    aligner.clustalExe = config.CLUSTAL
    aligner.msaExe = config.MSA


if len (sys.argv) != 3 or not os.path.isfile (sys.argv[2]):
    print 'usage: %s [species id] [fasta filepath]' % sys.argv[0]
    sys.exit(1)

init()

# The script takes a uniprot species ID and fasta filepath as input.
# Fasta headers are ignored, all sequences are assumed to be
# from the same species.
#
# A model is attempted to be built for every sequence in the fasta.

species = sys.argv [1]
seqs = parseFasta (open (sys.argv [2],'r'))

for ID in seqs:

    try:

        modeler.modelProc (seqs[ID], species)

    except:
        # Print the error to stdout, but don't quit.

        exc_type, exc_value, exc_traceback = sys.exc_info()
        stacktrace = ''.join(traceback.format_exception (exc_type, exc_value, exc_traceback))
        print stacktrace
