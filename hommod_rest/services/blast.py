from hommod_rest.services.modelutils import writeFasta, parseBlastXML
import logging
import os
import subprocess

from glob import glob

_log = logging.getLogger(__name__)


# Uses ncbi blast executables: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
class BlastService(object):

    def __init__(self):

        self.blastpExe = None
        self.templatesDB = None
        self.uniprotDB = None

    def _checkinit (self):

        from flask import current_app as flask_app

        if not self.templatesDB:
            self.templatesDB = flask_app.config ['TEMPLATESDB']

        if not self.uniprotDB:
            self.uniprotDB = flask_app.config ['UNIPROTDB']

        if not self.blastpExe:
            self.blastpExe = flask_app.config ['BLASTP']

    # Blast against a database of templates:
    def templateBlast (self, seq):

        self._checkinit ()
        _log.info ("performing template blast for\n%s" %seq)

        if not self.templatesDB:
            _log.error ("templates blast database not set")
            raise Exception("templates blast database not set")

        return self._blast (seq, self.templatesDB)

    # Blast against a proteome of a given species:
    def speciesBlast (self, seq, species):

        self._checkinit ()
        _log.info ("performing species blast for\n%s" %seq)

        if not self.uniprotDB:
            _log.error ("Species database directory not set")
            raise Exception("Species database directory not set")

        hits = self._blast(seq, self.uniprotDB)

        species_hits = {}
        for hitID in hits:
            if hitID.endswith ('_' + species):
                species_hits [hitID] = hits [hitID]

        return species_hits


    def _blast(self, querySeq, db):

        self._checkinit ()
        if not self.blastpExe:

            _log.error ('blastp executable not set')
            raise Exception ('blastp executable not set')

        queryFile = '/tmp/query%i.fasta' % (os.getpid())
        writeFasta({'query': querySeq}, queryFile)

        outFile = '/tmp/results%i.xml' % (os.getpid())

        # Run blastp and make sure it gives xml output format:
        cs = [self.blastpExe, '-query', queryFile, '-db', db,
              '-outfmt', '5', '-out', outFile]
        _log.debug(str(cs))

        try:
            feedback = subprocess.check_output (cs, stderr=subprocess.STDOUT)

            xmlstr = open(outFile, 'r').read ()

        except subprocess.CalledProcessError as e:

            _log.error ("blast failed:\n%s" % e.output)
            raise Exception ("blast failed:\n%s" % e.output)

        finally:
            if os.path.isfile (queryFile):
                os.remove (queryFile)

            if os.path.isfile (outFile):
                os.remove (outFile)

        if len (xmlstr) <= 0:

            _log.error ("no blast output:\n%s" % (feedback))
            raise Exception ("no blast output:\n%s" % (feedback))

        hits = parseBlastXML (xmlstr)

        _log.info ("found %d hits for %s in %s" %(len (hits), querySeq, db))

        return hits

blaster = BlastService()
