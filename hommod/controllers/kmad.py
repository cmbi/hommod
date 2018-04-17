import tempfile
import os
import subprocess
import logging

from hommod.models.align import TargetTemplateAlignment
from hommod.controllers.fasta import write_fasta, parse_fasta
from hommod.models.error import InitError


_log = logging.getLogger(__name__)

class KmadAligner:
    def __init__(self, kmad_exe=None):
        self.kmad_exe = kmad_exe

    def align(self, template_sequence, template_secstr, target_sequence,
              gap_open=-13.0, gap_extend=-0.4, modifier=3.0):

        _log.debug("kmad align {} {} {}".format(template_sequence, template_secstr, target_sequence))

        # Prevent kmad from adding insertions in bulges, replace those.
        template_secstr = self._remove_bulges(template_secstr, 'H', 3)
        template_secstr = self._remove_bulges(template_secstr, 'E', 3)

        if len(template_sequence) <= 0:
            raise ValueError("empty template sequence")
        if len(template_sequence) != len(template_secstr):
            raise ValueError("template sequence ({}) has different length than secondary structure ({})"
                             .format(len(template_sequence), len(template_secstr)))

        kmad_template_sequence = self._to_kmad_sequence(template_sequence, template_secstr)
        kmad_target_sequence = self._to_kmad_sequence(target_sequence)

        input_path = tempfile.mktemp()
        output_path = tempfile.mktemp()

        write_fasta(input_path, {'target': kmad_target_sequence,
                                 'template': kmad_template_sequence})
        try:
            self._run_kmad(input_path, output_path, gap_open, gap_extend, modifier)

            output_path += '_al'

            aligned = parse_fasta(output_path)
        finally:
            for path in [input_path, output_path]:
                if os.path.isfile(path):
                    os.remove(path)

        alignment = TargetTemplateAlignment(aligned['target'], aligned['template'])
        return alignment

    def _run_kmad(self, input_path, output_path, gap_open, gap_extend, modifier):

        if self.kmad_exe is None:
            raise InitError("kmad executable is not set")

        cmd = [self.kmad_exe, '-i', input_path, '-o', output_path,
               '-g', '%.1f' % gap_open, '-e', '%.1f' % gap_extend, '-s', '%.1f' % modifier, '-c']

        _log.debug(cmd)

        p = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.wait()

        if p.returncode != 0:
            raise RuntimeError(p.stderr.read())

    def _to_kmad_sequence(self, sequence, secstr=None):
        kmad_sequence = ''
        for i in range(len(sequence)):
            if secstr is not None and secstr[i] in ['H', 'E']:
                kmad_sequence += '%sA%sA' % (sequence[i], secstr[i])
            else:
                kmad_sequence += '%sAAA' % sequence[i]

        return kmad_sequence

    def _remove_bulges(self, secstr, type_, length):
        bulge_surrounding = type_ * length
        i = length
        while (i + length) < len(secstr):
            if secstr[i - length: i] == bulge_surrounding and \
                    secstr[i + 1: i + 1 + length] == bulge_surrounding:
                secstr = secstr[:i] + type_ + secstr[i + 1:]
                i += length
            i += 1

        return secstr

kmad_aligner = KmadAligner()
