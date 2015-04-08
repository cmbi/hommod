from hommod_rest.services.modelutils import writeFasta, parseBlastXML
import logging
import os
import subprocess

from glob import glob

_log = logging.getLogger(__name__)


class BlastService(object):

    def __init__(self):

        self.blastpExe = None
        self.templatesDB = None
        self.uniprotspeciesDir = None

    def templateBlast(self, seq):

        if not self.templatesDB:
            raise Exception("templates blast database not set")

        return _blast(seq, self.templatesDB)

    def speciesBlast(self, seq, species):

        if not self.uniprotspeciesDir:
            raise Exception("Species database directory not set")

        dbpath = os.path.join(self.uniprotspeciesDir, 'uniprot-%s' % species)
        if len(glob('%s.*' % dbpath)) <= 0:
                    raise Exception('Species database not found: '+dbpath)

        return _blast(seq, dbpath)


def _blast(querySeq, db):

    if not blaster.blastpExe:
        raise Exception('blastp executable not set')

    queryFile = '/tmp/query%i.fasta' % (os.getpid())
    writeFasta({'query': querySeq}, queryFile)

    outFile = '/tmp/results%i.xml' % (os.getpid())

    cs = [blaster.blastpExe, '-query', queryFile, '-db', db,
          '-outfmt', '5', '-out', outFile]
    _log.debug(str(cs))
    subprocess.call(cs)

    if not os.path.isfile(outFile):
        raise Exception("blast failed")

    hits = parseBlastXML(open(outFile, 'r').read())

    os.remove(queryFile)
    os.remove(outFile)

    return hits

blaster = BlastService()
