import os
from filelock import FileLock
from modelutils import idForSeq, writeFasta
import xml.etree.ElementTree as xmlElementTree
import bz2
import subprocess

import logging
_log = logging.getLogger(__name__)

class InterproDomain(object):

    # Represents a range of a sequence.

    def __init__(self, start, end, ac):

        self.start = start
        self.end = end
        self.ac = ac

# Uses interproscan to obtain data: https://code.google.com/p/interproscan/wiki/HowToDownload
class InterproService (object):

    def __init__(self):

        self.interproExe = None
        self.storageDir = None

    def _create_data_file(self, sequence):

        _log.info ("creating interpro file for sequence:\n%s" % sequence)

        if not self.interproExe or not self.storageDir:
            raise Exception("interproExe and storageDir must be set")

        if not os.path.isdir(self.storageDir):
            os.mkdir(self.storageDir)

        ID = idForSeq (sequence) # unique

        # We don't want two processes to build the same interpro file
        # at the same time. So use a lock file:
        lockfilepath = os.path.join(self.storageDir, 'lock_%s' % ID)

        lock = FileLock(lockfilepath)

        with lock:

            wd = os.getcwd()
            _in = os.path.join(wd, 'in%i.fasta' % os.getpid())
            out = os.path.join(self.storageDir, '%s.xml' % ID)

            # If the interpro file is already there, then don't
            # make it anew.
            if os.path.isfile(out + '.bz2'):
                return out + '.bz2'

            # Make the input file for interpro:
            writeFasta({ID: sequence}, _in)

            # Run interproscan and wait for it to finish:
            cmd = '%s --goterms --formats xml --input %s --outfile %s --seqtype p' % \
                  (self.interproExe, _in, out)
            
            p = subprocess.Popen (cmd, shell=True, stderr=subprocess.PIPE)
            p.wait ()

            os.remove(_in)

            if not os.path.isfile (out):

                # Get the console output:
                errstr = p.stderr.read ()

                _log.error ("interpro did not create %s:\n%s" %
                            (out, errstr))

                raise Exception('interprofile %s not created:\n%s' % (out, errstr))

            # hommod assumes the interpro files to be bzip2 compressed by default:
            os.system ('bzip2 -f %s' % out)
            out = out + '.bz2'

            _log.info ("interprofile %s sucessfully created" % out)

            return out


    # Creates an interpro file for the given sequence and parses it.
    def getInterproDomainLocations(self, sequence):

        _log.info ("getting interpro domains for sequence:\n%s" % sequence)

        filepath = self._create_data_file(sequence)
        ID = idForSeq(sequence)

        # Get the xml tag of type protein and the given sequence ID:
        root = xmlElementTree.fromstring (bz2.BZ2File (filepath, 'r').read())
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

        # Look for pattern matches in the file:
        matches = []
        for child in protein:
            if child.tag.endswith('}matches'):
                matches = child
            elif child.tag == 'match':
                matches.append(child)

        # Get the ranges of the matched parts of the sequence:
        domains = []
        for match in matches:

            # Look at other aspects of the match as well.
            # We have conditions, some matches aren't added.
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

                # If the range is very short, it might not be an actual domain.
                # So 'bShortDomain' must be set as an additional confirmation.
                if length > 20 or bShortDomain:
                    domains.append(InterproDomain(start, end, ac))

        _log.info ("successfully retrieved interpro domains for sequence:\n%s" % sequence)

        return domains

interpro = InterproService()
