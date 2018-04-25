import requests
import logging

from hommod.models.error import InitError, ServiceError
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

        r = requests.get(fasta_url)
        while r.status_code == 500:
            r = requests.get(fasta_url)

        r.raise_for_status()

        fa = parse_fasta_from_string(r.text)

        return fa.values()[0]


uniprot = UniprotService()
