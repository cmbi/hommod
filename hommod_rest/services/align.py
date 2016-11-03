import os
import subprocess

from hommod_rest.services.modelutils import (
    writeFasta, parseFasta, removeBulges)

import logging
_log = logging.getLogger(__name__)


class AlignService(object):

    def __init__(self, clustal_exe=None, kmad_exe=None):
        self._clustal_exe = clustal_exe
        self._kmad_exe = kmad_exe

    @property
    def clustal_exe(self):
        return self._clustal_exe

    @clustal_exe.setter
    def clustal_exe(self, clustal_exe):
        self._clustal_exe = clustal_exe

    @property
    def kmad_exe(self):
        return self._kmad_exe

    @kmad_exe.setter
    def kmad_exe(self, kmad_exe):
        self._kmad_exe = kmad_exe

    def _checkinit(self):
        if not self._clustal_exe:
            raise Exception("clustal_exe not set")
        if not self._kmad_exe:
            raise Exception("kmad_exe not set")

    def clustal_align(self, d):
        """
        Uses clustalw2 for alignment: http://www.clustal.org/clustal2/
        Simply uses dictionaries for input and output.
        """

        self._checkinit()

        _log.info("clustal aligning %s" % str(d))

        wd = os.path.abspath(os.getcwd())

        for key in d.keys():
            if '|' in key:
                _log.error('Invalid syntax for key: ' + key)
                raise Exception('Invalid syntax for key: ' + key)
            if len(d[key]) == 0:
                _log.error('empty sequence for %s' % key)
                raise Exception('empty sequence for %s' % key)

        _in = os.path.join(wd, 'in%i.fasta' % os.getpid())
        writeFasta(d, _in)

        out = os.path.join(wd, 'out%i.fasta' % os.getpid())

        args = [
            self.clustal_exe, '-TYPE=PROTEIN', '-OUTPUT=FASTA',
            '-PWMATRIX=BLOSUM', '-OUTFILE=%s' % out, '-INFILE=%s' % _in
        ]

        try:
            subprocess.call(args, stdout=subprocess.PIPE)

        finally:
            DND = os.path.splitext(_in)[0] + '.dnd'
            if os.path.isfile(DND):
                os.remove(DND)
            if os.path.isfile(_in):
                os.remove(_in)

        try:
            d = parseFasta(open(out, 'r'))
        finally:
            if os.path.isfile(out):
                os.remove(out)

        _log.info("successfully created a clustal alignment")

        return d

    def kmad_align(self, pdbSeq, pdbSecStr, tarSeq,
                  gapOpen=-13.0, gapExt=-0.4, modifier=3.0):
        """
        Uses Joanna Lange's alignment program, requires secondary structure
        information. Output is a dictionary with: 'template' and 'target' as
        ids, pointing to sequences.
        """

        self._checkinit()
        _log.info("making pairwise kmad alignment")

        # Prevent kmad from adding insertions in bulges.
        # Fool the program by making
        # it think they're helix/strand residues:
        pdbSecStr = removeBulges(pdbSecStr, 'H', 3)
        pdbSecStr = removeBulges(pdbSecStr, 'E', 3)

        if len(pdbSeq) == 0:
            _log.error('empty pdb seq')
            raise Exception('empty pdb seq')

        if len(pdbSeq) != len(pdbSecStr):
            _log.error('pdb seq and pdb secstr are not of same length')
            raise Exception('pdb seq and pdb secstr are not of same length')

        seq1 = ''
        for i in range(len(pdbSeq)):
            aa = pdbSeq[i]
            ss = pdbSecStr[i]
            if ss in ['H', 'E']:
                seq1 += '%sA%sA' % (aa, ss)
            else:
                seq1 += '%sAAA' % aa

        seq2 = ''
        for aa in tarSeq:
            seq2 += '%sAAA' % aa

        toalignpath = 'toalign%i.fasta' % os.getpid()
        alignedpath = 'aligned%i' % os.getpid()

        _input = '>template\n%s\n>target\n%s\n' % (seq1, seq2)
        open(toalignpath, 'w').write(_input)

        cmd = '%s -i %s -o %s -g %.1f -e %.1f -s %.1f -c; exit 0' % \
              (self.kmad_exe, toalignpath, alignedpath,
               gapOpen, gapExt, modifier)

        try:
            feedback = subprocess.check_output(cmd, shell=True,
                                               stderr=subprocess.STDOUT)

        finally:
            if os.path.isfile(toalignpath):
                os.remove(toalignpath)

        alignedpath += '_al'

        if not os.path.isfile(alignedpath):

            error = ('alignment file %s not created,' + 
                     'kmad input:\n%s\nkmad error:\n%s') % \
                            (alignedpath, _input, feedback)
            _log.error(error)
            raise Exception(error)

        try:
            aligned = parseFasta(open(alignedpath, 'r'))

        finally:
            if os.path.isfile(alignedpath):
                os.remove(alignedpath)

        if aligned['template'].replace('-', '') != pdbSeq:
            error = 'kmad output mismatch:\npdbSeq:' + pdbSeq + \
                            'aligned:' + aligned['template']
            _log.error(error)
            raise Exception(error)

        _log.debug("successfully created a pairwise kmad alignment")

        return aligned

aligner = AlignService()
