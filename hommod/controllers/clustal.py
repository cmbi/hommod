import tempfile
import os
import subprocess

from hommod.models.align import Alignment
from hommod.controllers.fasta import parse_fasta, write_fasta
from hommod.models.error import InitError


class ClustalAligner:

    def __init__(self, clustalw_exe=None):
        self.clustalw_exe = clustalw_exe

    def align(self, input_):
        if self.clustalw_exe is None:
            raise InitError("clustalw executable is not set")

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

            return Alignment(parse_fasta(output_path))
        finally:
            for path in [input_path, output_path]:
                if os.path.isfile(path):
                    os.remove(path)


clustal_aligner = ClustalAligner()
