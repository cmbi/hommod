import logging
import datetime
from time import time, sleep
import xml.etree.ElementTree as ET
import requests

from hommod.services.helpers.cache import cache_manager as cm
from hommod.models.range import SequenceRange
from hommod.models.error import InitError, ServiceError


_log = logging.getLogger(__name__)

class InterproService:
    def __init__(self, url=None, email=None,
                 job_timeout=60*60*2, http_timeout=60, poll_interval=10):
        self.url = url
        self.email = email
        self.http_timeout = http_timeout
        self.job_timout = job_timeout
        self.poll_interval = poll_interval

    @cm.cache()
    def get_domain_ranges(self, sequence):
        if self.url is None:
            raise InitError("interpro url is not set")

        job_id = self._interpro_submit(sequence)

        t0 = time()
        while (time() - t0) < self.job_timout:

            status = self._interpro_status(job_id)

            if status in ['RUNNING', 'PENDING', 'STARTED', 'QUEUED']:
                sleep(self.poll_interval)
            elif status == 'NOT_FOUND':
                job_id = self._interpro_submit(sequence)
            else:
                break

        if status == 'RUNNING':
            raise ServiceError("inteproscan job timed out")
        elif status in ['FAILURE', 'ERROR']:

            response_text = self._interpro_result(job_id)
            if InterproService._is_usable_output(response_text):
                return self._parse_interpro_ranges(response_text)

            raise ServiceError(self._interpro_error(job_id))
        elif status != 'FINISHED':
            raise ServiceError("inteproscan job status = " + status)

        xml_str = self._interpro_result(job_id)

        return self._parse_interpro_ranges(xml_str)


    def _interpro_submit(self, sequence):

        params = {'email': self.email,
                  'sequence': sequence,
                  'goterms': True,
                  'pathways': False}

        submit_url = '/'.join([self.url, 'run'])
        try:
            r = requests.post(submit_url, data=params, timeout=self.http_timeout)
        except requests.exceptions.ConnectTimeout:
            raise ServiceError("timeout connecting with interpro")

        r.raise_for_status()

        return r.text

    def _interpro_status(self, job_id):

        status_url = '/'.join([self.url, 'status', job_id])
        try:
            r = requests.get(status_url, timeout=self.http_timeout)
        except requests.exceptions.ConnectTimeout:
            raise ServiceError("timeout connecting with interpro")

        r.raise_for_status()

        return r.text

    def _interpro_result(self, job_id):

        result_url = '/'.join([self.url, 'result', job_id, 'xml'])
        try:
            r = requests.get(result_url, timeout=self.http_timeout)
        except requests.exceptions.ConnectTimeout:
            raise ServiceError("timeout connecting with interpro")

        r.raise_for_status()

        return r.text

    def _interpro_error(self, job_id):

        result_url = '/'.join([self.url, 'result', job_id, 'xml'])
        try:
            r = requests.get(result_url, timeout=self.http_timeout)
        except requests.exceptions.ConnectTimeout:
            raise ServiceError("timeout connecting with interpro")

        r.raise_for_status()

        return r.text

    def _parse_interpro_ranges(self, xml_str):

        ranges = []

        ns_map = {}
        if "/software/unix/" in xml_str:
            ns_map["p"] = "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/schemas"
        else:
            ns_map["p"] = "http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5"

        root = ET.fromstring(xml_str)
        sequence_elem = root.find('./p:protein/p:sequence', namespaces=ns_map)
        sequence = sequence_elem.text

        for match_elem in root.find('./p:protein/p:matches', namespaces=ns_map):
            entry_elem = match_elem.find('.//p:signature/p:entry', namespaces=ns_map)
            if entry_elem is None:
                continue

            ac = entry_elem.get('ac')
            desc = entry_elem.get('desc')

            allow_short_domain = self._short_domain_allowed(entry_elem)

            for location_elem in match_elem.find('.//p:locations', namespaces=ns_map):
                start = int(location_elem.get('start')) - 1
                end = int(location_elem.get('end')) - 1
                length = end - start

                if length > 20 or allow_short_domain:
                    range_ = SequenceRange(start, end, sequence)
                    range_.ac = ac
                    ranges.append(range_)

        return ranges

    def _short_domain_allowed(self, entry_elem):

        if 'zinc finger' in entry_elem.get('desc').lower():
            return True

        return False

    @staticmethod
    def _is_usable_output(response_text):
        try:
            x = ET.fromstring(response_text)

        except ET.ParseError:
            return False

        _log.debug("x.tag is {}".format(x.tag))

        return x.tag == "{http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5}protein-matches"


interpro = InterproService()
