import os
from filelock import FileLock
from modelutils import idForSeq, writeFasta
import xml.etree.ElementTree as xmlElementTree
import bz2


class InterproDomain(object):

        def __init__(self, start, end, ac):

                self.start = start
                self.end = end
                self.ac = ac

# Uses interproscan to obtain data: https://code.google.com/p/interproscan/wiki/HowToDownload
class InterproService(object):

    def __init__(self):

        self.interproExe = None
        self.storageDir = None

    def _create_data_file(self, sequence):

        if not self.interproExe or not self.storageDir:
            raise Exception("interproExe and storageDir must be set")

        if not os.path.isdir(self.storageDir):
            os.mkdir(self.storageDir)

        ID = idForSeq(sequence)

        lockfilepath = os.path.join(self.storageDir, 'lock_%s' % ID)

        lock = FileLock(lockfilepath)

        with lock:

            wd = os.getcwd()
            _in = os.path.join(wd, 'in%i.fasta' % os.getpid())
            out = os.path.join(self.storageDir, '%s.xml' % ID)

            if os.path.isfile(out + '.bz2'):
                return out + '.bz2'

            writeFasta({ID: sequence}, _in)

            os.system(('%s --goterms --formats xml ' +
                       '--input %s --outfile %s --seqtype p')
                      % (self.interproExe, _in, out))

            os.remove(_in)

            if not os.path.isfile(out):

                raise Exception('interprofile not created: %s' % out)

            os.system('bzip2 -f %s' % out)
            return out + '.bz2'


    # Creates data file if it doesn't exist yet.
    def getInterproDomainLocations(self, sequence):

        filepath = self._create_data_file(sequence)
        ID = idForSeq(sequence)

        root = xmlElementTree.fromstring(bz2.BZ2File(filepath, 'r').read())
        for child in root:
            if child.tag.endswith('}protein'):
                p = child
                for child in p:
                    if child.tag.endswith('}xref') and \
                            child.attrib['id'] == ID:
                        protein = p
                        break
            elif child.tag == 'protein' and child.attrib['id'] == ID:
                protein = child

        if protein is None:
            raise Exception('Protein not found: '+ID)

        matches = []
        for child in protein:
            if child.tag.endswith('}matches'):
                matches = child
            elif child.tag == 'match':
                matches.append(child)

        domains = []
        for match in matches:
            bShortDomain = False
            ac = None
            locations = []
            for child in match:
                if child.tag.endswith('}signature'):
                    for c in child:
                        if c.tag.endswith('}entry'):
                            desc = c.attrib['desc'].lower()
                            ac = c.attrib['ac']

                            # only sinc finger domains are allowed to be short
                            if 'zinc finger' in desc:
                                bShortDomain = True

                if child.tag == 'lcn':
                    locations.append(child)
                if child.tag == 'ipr' and \
                        'zinc finger' in child.attrib['name'].lower():
                    bShortDomain = True

                if child.tag.endswith('}locations'):
                    locations = child

            for location in locations:
                start = int(location.attrib['start'])-1
                end = int(location.attrib['end'])-1
                length = end - start

                if length > 20 or bShortDomain:
                    domains.append(InterproDomain(start, end, ac))

        return domains

interpro = InterproService()
