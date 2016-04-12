import os
import subprocess

from hommod_rest.services.modelutils import (
    writeFasta, parseFasta, removeBulges)

import logging
_log = logging.getLogger(__name__)

class alignService(object):

    def __init__(self):

        self.clustalExe = None
        self.msaExe = None

    def _checkinit (self):

        from flask import current_app as flask_app

        self.clustalExe = flask_app.config ['CLUSTAL']
        self.msaExe = flask_app.config ['MSA']

    # Use Clustalw2 for alignment: http://www.clustal.org/clustal2/
    # Simply uses dictionaries for input and output
    def clustalAlign(self, d):

        self._checkinit ()

        _log.info ("clustal aligning %s" % str (d))

        if not self.clustalExe:
            _log.error ("clustal path not set")
            raise Exception ("clustal path not set")

        wd = os.path.abspath(os.getcwd())

        for key in d.keys():
            if '|' in key:
                _log.error ('Unallowed syntax for key: ' + key)
                raise Exception('Unallowed syntax for key: ' + key)
            if len(d[key]) == 0:
                _log.error ('empty sequence for %s' % key)
                raise Exception('empty sequence for %s' % key)

        _in = os.path.join(wd, 'in%i.fasta' % os.getpid())
        writeFasta(d, _in)

        out = os.path.join(wd, 'out%i.fasta' % os.getpid())

        args = [
            self.clustalExe, '-TYPE=PROTEIN', '-OUTPUT=FASTA',
            '-PWMATRIX=BLOSUM', '-OUTFILE=%s' % out, '-INFILE=%s' % _in
        ]
        print args

        try:
            subprocess.call(args, stdout=subprocess.PIPE)

        finally:
            os.remove(os.path.splitext(_in)[0] + '.dnd')
            os.remove(_in)

        try:
            d = parseFasta(open(out, 'r'))
        finally:
            os.remove(out)

        _log.info ("successfully created a clustal alignment")

        return d

    # Use Joanna Lange's alignment program, requires secondary structure information.
    # Output is a dictionary with: 'template' and 'target' as ids, pointing to sequences.
    def msaAlign(self, pdbSeq, pdbSecStr, tarSeq,
                 gapOpen=-13.0, gapExt=-0.4, modifier=3.0):

        self._checkinit ()
        _log.info ("making pairwise MSA alignment")

        if not self.msaExe:
            raise Exception("msa path not set")

        # Prevent MSA from adding insertions in bulges.
        # Fool the program by making
        # it think they're helix/strand residues:
        pdbSecStr = removeBulges (pdbSecStr, 'H', 3)
        pdbSecStr = removeBulges (pdbSecStr, 'E', 3)

        if len(pdbSeq) == 0:
            _log.error ('empty pdb seq')
            raise Exception('empty pdb seq')

        if len(pdbSeq) != len(pdbSecStr):
            _log.error ('pdb seq and pdb secstr are not of same length')
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
        open(toalignpath, 'w').write('>template\n%s\n>target\n%s\n' %
                                     (seq1, seq2))

        cmd = '%s -i %s -o %s -g %.1f -e %.1f -s %.1f -c' % \
              (self.msaExe, toalignpath, alignedpath,
               gapOpen, gapExt, modifier)

        p = subprocess.Popen (cmd, shell=True, stderr=subprocess.PIPE)
        p.wait ()

        alignedpath += '_al'

        if not os.path.isfile (alignedpath):

            errstr = p.stderr.read ()

            _log.error ('alignment file %s not created, MSA error:\n%s' %
                        (alignedpath, errstr))
            raise Exception ('alignment file %s not created, MSA error:\n%s' %
                             (alignedpath, errstr))

        aligned = parseFasta (open(alignedpath, 'r'))

        os.remove (alignedpath)
        os.remove (toalignpath)

        if aligned['template'].replace('-', '') != pdbSeq:
            _log.error ('MSA output mismatch:\npdbSeq:' + pdbSeq +
                            'aligned:' + aligned['template'])
            raise Exception('MSA output mismatch:\npdbSeq:' + pdbSeq +
                            'aligned:' + aligned['template'])

        _log.info ("successfully created a pairwise MSA alignment")

        return aligned

aligner = alignService()
