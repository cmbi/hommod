import tempfile
import os
import subprocess

from hommod.models.align import Alignment
from hommod.controllers.fasta import parse_fasta, write_fasta
from hommod.models.error import InitError


class ClustalAligner:

    def __init__(self, clustalw_exe=None):
        self.clustalw_exe = clustalw_exe

    def _lowercase_escape(self, key):
        "makes lowercase distinctive from uppercase, for clustal"

        s = ""
        for i in range(len(key)):
            if key[i].islower() or key[i] == "_":
                s += '_' + key[i]
            else:
                s += key[i]

        return s

    def _lowercase_unescape(self, esckey):
        "undoes _lowercase_escape"

        s = ""
        i = 0
        while i < len(esckey):
            if esckey[i] == '_':
                s += esckey[i + 1].lower()
                i += 2
            else:
                s += esckey[i]
                i += 1

        return s

    def _fix_input(self, input_):
        "makes sure clustal can work with it"

        inp = {}
        for key in input_:
            inp[self._lowercase_escape(key)] = input_[key]

        return inp

    def _fix_output(self, output):
        "undoes _fix_input"

        outp = {}
        for esckey in output:
            outp[self._lowercase_unescape(esckey)] = output[esckey]
        return outp

    def align(self, input_):
        if self.clustalw_exe is None:
            raise InitError("clustalw executable is not set")

        input_ = self._fix_input(input_)

        input_path = tempfile.mktemp()
        output_path = tempfile.mktemp()

        write_fasta(input_path, input_)

        cmd = [self.clustalw_exe, '-TYPE=PROTEIN', '-OUTPUT=FASTA',
               '-PWMATRIX=BLOSUM', '-OUTFILE=%s' % output_path, '-INFILE=%s' % input_path]

        try:
            p = subprocess.Popen(cmd,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p.wait()

            if p.returncode != 0:
                raise RuntimeError("%s for %s" % (p.stderr.read().decode('ascii'), str(input_)))

            return Alignment(self._fix_output(parse_fasta(output_path)))
        finally:
            for path in [input_path, output_path]:
                if os.path.isfile(path):
                    os.remove(path)


clustal_aligner = ClustalAligner()
