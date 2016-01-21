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
        self.uniprotspeciesDir = None

    # Blast against a database of templates:
    def templateBlast (self, seq):

        _log.info ("performing template blast for\n%s" %seq)

        if not self.templatesDB:
            _log.error ("templates blast database not set")
            raise Exception("templates blast database not set")

        return _blast (seq, self.templatesDB)

    # Blast against a proteome of a given species:
    def speciesBlast (self, seq, species):

        _log.info ("performing species blast for\n%s" %seq)

        if not self.uniprotspeciesDir:
            _log.error ("Species database directory not set")
            raise Exception("Species database directory not set")

        dbpath = os.path.join(self.uniprotspeciesDir, 'uniprot-%s' % species.upper ())
        if len(glob('%s.*' % dbpath)) <= 0:
                    _log.error ('Species database not found: ' + dbpath)
                    raise Exception('Species database not found: ' + dbpath)

        return _blast(seq, dbpath)


def _blast(querySeq, db):

    if not blaster.blastpExe:

        _log.error ('blastp executable not set')
        raise Exception('blastp executable not set')

    queryFile = '/tmp/query%i.fasta' % (os.getpid())
    writeFasta({'query': querySeq}, queryFile)

    outFile = '/tmp/results%i.xml' % (os.getpid())

    # Runt blastp and make sure it gives xml output format:
    cs = [blaster.blastpExe, '-query', queryFile, '-db', db,
          '-outfmt', '5', '-out', outFile]
    _log.debug(str(cs))

    p = subprocess.Popen (cs, stderr=subprocess.PIPE)
    p.wait ()

    if not os.path.isfile (outFile):

        errstr = p.stderr.read ()

        _log.error ("blast failed:\n%s" % errstr)
        raise Exception ("blast failed:\n%s" % errstr)

    hits = parseBlastXML (open(outFile, 'r').read())

    os.remove (queryFile)
    os.remove (outFile)

    _log.info ("found %d hits for %s in %s" %(len (hits), querySeq, db))

    return hits

blaster = BlastService()
