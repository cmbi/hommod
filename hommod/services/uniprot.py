from urllib2 import urlopen
import logging

from hommod.models.error import InitError
from hommod.controllers.fasta import parse_fasta_from_string


_log = logging.getLogger(__name__)

class UniprotService:
    def __init__(self, url=None):
        self.url = url

    def get_sequence(self, ac):
        if self.url is None:
            raise InitError("uniprot url is not set")

        fasta_url = self.url + '/' + ac + '.fasta'

        _log.debug(fasta_url)

        fa = parse_fasta_from_string(urlopen(fasta_url).read())

        return fa.values()[0]


uniprot = UniprotService()
