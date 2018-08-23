import os
import logging
import tempfile
import subprocess
import xml.etree.ElementTree as ET

from hommod.models.error import InitError, RecoverableError
from hommod.controllers.fasta import write_fasta
from hommod.models.align import BlastAlignment

_log = logging.getLogger(__name__)


class Blaster:

    def __init__(self, blastp_exe=None):
        self.blastp_exe = blastp_exe

    def blastp(self, sequence, databank):
        if self.blastp_exe is None:
            raise InitError("blastp executable is not set")

        input_path = tempfile.mktemp()
        output_path = tempfile.mktemp()

        write_fasta(input_path, {'target': sequence})

        cmd = [self.blastp_exe, '-query', input_path, '-db', databank,
               '-outfmt', '5', '-out', output_path]

        _log.debug("{}".format(cmd))

        try:
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p.wait()

            if p.returncode != 0:
                err_msg = p.stderr.read()
                if err_msg.startswith("BLAST Database error: No alias or index file found for protein database"):
                    raise RecoverableError(err_msg)

                raise RuntimeError("%s for databank %s, sequence %s"
                                   % (err_msg, databank, sequence))

            with open(output_path, 'r') as f:
                xml_str = f.read()
        finally:
            for path in [input_path, output_path]:
                if os.path.isfile(path):
                    os.remove(path)

        return self._parse_alignments(xml_str, sequence, databank)

    def _parse_alignments(self, xml_str, full_query_sequence, databank):
        hits = {}
        root = ET.fromstring(xml_str)
        iterations = root.find('BlastOutput_iterations')
        for it in iterations.findall('Iteration'):
            for mem in it.findall('Iteration_hits'):
                for hit in mem.findall('Hit'):
                    hit_id = hit.find('Hit_def').text
                    hits[hit_id] = []

                    hsps = hit.find('Hit_hsps')
                    for hsp in hsps.findall('Hsp'):

                        query_start = int(hsp.find('Hsp_query-from').text)
                        query_end = int(hsp.find('Hsp_query-to').text)
                        query_alignment = hsp.find('Hsp_qseq').text

                        subject_start = int(hsp.find('Hsp_hit-from').text)
                        subject_end = int(hsp.find('Hsp_hit-to').text)
                        subject_alignment = hsp.find('Hsp_hseq').text

                        hits[hit_id].append(BlastAlignment(hit_id,
                                                           full_query_sequence,
                                                           databank,
                                                           query_start,
                                                           query_end,
                                                           query_alignment,
                                                           subject_start,
                                                           subject_end,
                                                           subject_alignment))
        return hits

blaster = Blaster()
