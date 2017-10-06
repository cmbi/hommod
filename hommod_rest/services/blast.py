from hommod_rest.services.modelutils import writeFasta, parseBlastXML
import logging
import os
import subprocess

_log = logging.getLogger(__name__)


# Uses ncbi blast executables: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
class BlastService(object):
    def __init__(self, blastp_exe=None, templates_db=None, uniprot_db=None, models_db=None):
        self._blastp_exe = blastp_exe
        self._templates_db = templates_db
        self._uniprot_db = uniprot_db
        self._models_db = models_db

    @property
    def blastp_exe(self):
        return self._blastp_exe

    @blastp_exe.setter
    def blastp_exe(self, blastp_exe):
        self._blastp_exe = blastp_exe

    @property
    def templates_db(self):
        return self._templates_db

    @templates_db.setter
    def templates_db(self, templates_db):
        self._templates_db = templates_db

    @property
    def uniprot_db(self):
        return self._uniprot_db

    @uniprot_db.setter
    def uniprot_db(self, uniprot_db):
        self._uniprot_db = uniprot_db

    @property
    def models_db(self):
        return self._models_db

    @models_db.setter
    def models_db(self, models_db):
        self._models_db = models_db

    def _checkinit(self):
        if not self.templates_db:
            raise Exception("templates_db not set")

        if not self.uniprot_db:
            raise Exception("uniprot_db not set")

        if not self.blastp_exe:
            raise Exception("blastp_exe not set")

        if not self.models_db:
            raise Exception("models_db not set")

    def blast_models(self, seq):
        """
        Blast against a database of models
        """

        self._checkinit()
        _log.info("performing model blast for\n%s" %seq)

        return self._blast(seq, self.models_db)

    def blast_templates(self, seq):
        """
        Blast against a database of templates.
        """

        self._checkinit()
        _log.info("performing template blast for\n%s" %seq)

        return self._blast(seq, self.templates_db)

    def blast_species(self, seq, species):
        """
        Blast against a proteome of a given species.
        """

        self._checkinit()
        _log.info("performing species blast for\n%s" %seq)

        hits = self._blast(seq, self.uniprot_db)

        species_hits = {}
        for hitID in hits:
            if hitID.endswith('_' + species):
                species_hits[hitID] = hits[hitID]

        return species_hits

    def _blast(self, querySeq, db):
        self._checkinit()

        queryFile = '/tmp/query%i.fasta' % (os.getpid())
        writeFasta({'query': querySeq}, queryFile)

        outFile = '/tmp/results%i.xml' % (os.getpid())

        # Run blastp and make sure it gives xml output format:
        cs = [self.blastp_exe, '-query', queryFile, '-db', db,
              '-outfmt', '5', '-out', outFile]
        _log.debug(str(cs))

        try:
            feedback = subprocess.check_output(cs, stderr=subprocess.STDOUT)
            xmlstr = open(outFile, 'r').read()
        except subprocess.CalledProcessError as e:
            _log.error("blast failed:\n%s" % e.output)
            raise Exception("blast failed:\n%s" % e.output)
        finally:
            if os.path.isfile(queryFile):
                os.remove(queryFile)

            if os.path.isfile(outFile):
                os.remove(outFile)

        if len(xmlstr) <= 0:
            _log.error("no blast output:\n%s" % (feedback))
            raise Exception("no blast output:\n%s" % (feedback))

        hits = parseBlastXML(xmlstr)

        _log.info("found %d hits for %s in %s" %(len(hits), querySeq, db))

        return hits

blaster = BlastService()
