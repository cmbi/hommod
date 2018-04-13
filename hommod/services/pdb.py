import logging
from urllib2 import Request, urlopen
from io import BytesIO
from gzip import GzipFile


_log = logging.getLogger(__name__)

def get_pdb_contents(pdbid):
    part = pdbid[1:3].lower()
    pdb_url = (
        'ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/%s/pdb%s.ent.gz'
        % (part, pdbid.lower()))

    _log.debug(pdb_url)

    request = Request(pdb_url)
    request.add_header('Accept-encoding', 'gzip')

    response = urlopen(request)

    buf = BytesIO(response.read())
    return GzipFile(fileobj=buf).read()
