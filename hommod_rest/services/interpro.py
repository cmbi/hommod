import os
from filelock import FileLock
from modelutils import idForSeq
import xml.etree.ElementTree as xmlElementTree
import bz2
import urllib
import urllib2
import platform
from time import sleep
from glob import glob

import logging
_log = logging.getLogger(__name__)

_interpro_base_url = 'http://www.ebi.ac.uk/Tools/services/rest/iprscan5'


def _interpro_user_agent():

    # Agent string for urllib2 library.
    urllib_agent = 'Python-urllib/%s' % urllib2.__version__

    clientRevision = '$Revision: 2809 $'
    clientVersion = '0'
    if len(clientRevision) > 11:
        clientVersion = clientRevision[11:-2]

    # Prepend client specific agent string.
    user_agent = 'EBI-Sample-Client/%s (%s; Python %s; %s) %s' % (

        clientVersion, os.path.basename(__file__),
        platform.python_version(), platform.system(),
        urllib_agent
    )

    return user_agent


def _interpro_post(url, data):

    requestData = urllib.urlencode(data)

    # Set the HTTP User-agent.
    user_agent = _interpro_user_agent()
    http_headers = {'User-Agent': user_agent}
    req = urllib2.Request(url, None, http_headers)

    # Make the submission (HTTP POST).
    reqH = urllib2.urlopen(req, requestData)
    result = reqH.read()
    reqH.close()

    return result


def _interpro_get(url):

    user_agent = _interpro_user_agent()
    http_headers = {'User-Agent': user_agent}
    req = urllib2.Request(url, None, http_headers)

    # Make the request (HTTP GET).
    reqH = urllib2.urlopen(req)
    result = reqH.read()
    reqH.close()

    return result


def _interpro_submit(sequence):

    params = {
        'email': "c.baakman@radboudumc.nl",
        'sequence': sequence,
        'goterms': True,
        'pathways': False
    }

    return _interpro_post(_interpro_base_url + '/run/', params)


def _interpro_get_status(jobid):
    return _interpro_get(_interpro_base_url + '/status/' + jobid)


def _interpro_get_result(jobid):
    return _interpro_get(_interpro_base_url + '/result/' + jobid + '/xml')


class InterproDomain(object):

    # Represents a range of a sequence.

    def __init__(self, start, end, ac):

        self.start = start
        self.end = end
        self.ac = ac

_MAX_INTERPRO_JOBS = 30


# Uses interproscan to obtain data: https://code.google.com/p/interproscan/wiki/HowToDownload
class InterproService (object):
    def __init__(self, storage_dir=None):
        self._storage_dir = storage_dir

    @property
    def storage_dir(self):
        return self._storage_dir

    @storage_dir.setter
    def storage_dir(self, storage_dir):
        self._storage_dir = storage_dir

    def _checkinit(self):
        if not self._storage_dir:
            raise Exception("storage_dir not set")

    def _interpro_lockfile_path(self, ID):

        return os.path.join(self.storage_dir, 'lock_%s' % ID)

    def _interpro_count_lockfiles(self):

        return len(glob(os.path.join(self.storage_dir, "lock_*")))

    def _create_data_file(self, sequence):

        self._checkinit()
        _log.info("creating interpro file for sequence:\n%s" % sequence)

        if not os.path.isdir(self.storage_dir):
            os.mkdir(self.storage_dir)

        ID = idForSeq(sequence)  # unique

        outfilepath = os.path.join(self.storage_dir, '%s.xml.bz2' % ID)
        if os.path.isfile(outfilepath):
            return outfilepath

        # We don't want two processes to build the same interpro file
        # at the same time. So use a lock file:
        lock = FileLock(self._interpro_lockfile_path(ID))

        with lock:
            # Wait for a place in line to start an interpro job
            while self._interpro_count_lockfiles() >= _MAX_INTERPRO_JOBS:
                sleep(10)

            # Wait for the interpro server to finish:
            jobid = _interpro_submit(sequence)
            while True:

                status = _interpro_get_status(jobid)
                _log.debug("intepro job status: " + status)

                if status in ['RUNNING', 'PENDING', 'STARTED']:

                    sleep(10)
                else:
                    break

            if status != 'FINISHED':
                raise Exception("inteproscan job status = " + status)

            xmlstr = _interpro_get_result(jobid)

            # Write results to file
            bz2.BZ2File(outfilepath, 'w').write(xmlstr)

            return outfilepath

    # Creates an interpro file for the given sequence and parses it.
    def getInterproDomainLocations(self, sequence):

        self._checkinit()
        _log.info("getting interpro domains for sequence:\n%s" % sequence)

        filepath = self._create_data_file(sequence)
        ID = idForSeq(sequence)

        # Get the xml tag of type protein and the given sequence ID:
        protein = None
        s = ''
        root = xmlElementTree.fromstring(bz2.BZ2File(filepath, 'r').read())
        for p in root.iter():
            if p.tag.endswith('}protein') or p.tag == 'protein':

                protein = p
                break

        if protein is None:
            raise Exception('Protein not found: ' + ID + ' in' + s)

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

                start = int(location.attrib['start']) - 1
                end = int(location.attrib['end']) - 1
                length = end - start

                # If the range is very short, it might not be an actual domain.
                # So 'bShortDomain' must be set as an additional confirmation.
                if length > 20 or bShortDomain:
                    domains.append(InterproDomain(start, end, ac))

        _log.info("successfully retrieved interpro domains for sequence:\n%s" % sequence)

        return domains

interpro = InterproService()
