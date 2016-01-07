import os
import subprocess

from hommod_rest.services.modelutils import (
    writeFasta, parseFasta, removeBulges)


class alignService(object):

    def __init__(self):

        self.clustalExe = None
        self.msaExe = None

    # Use Clustalw2 for alignment: http://www.clustal.org/clustal2/
    # Simply uses dictionaries for input and output
    def clustalAlign(self, d):

        if not self.clustalExe:
            raise Exception("clustal path not set")

        wd = os.path.abspath(os.getcwd())

        for key in d.keys():
            if '|' in key:
                raise Exception('Unallowed syntax for key: ' + key)
            if len(d[key]) == 0:
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

        return d

    # Use Joanna Lange's alignment program, requires secondary structure information.
    # Output is a dictionary with: 'template' and 'target' as ids, pointing to sequences.
    def msaAlign(self, pdbSeq, pdbSecStr, tarSeq,
                 gapOpen=-13.0, gapExt=-0.4, modifier=3.0):

        if not self.msaExe:
            raise Exception("msa path not set")

        # Prevent MSA from adding insertions in bulges.
        # Fool the program by making
        # it think they're helix/strand residues:
        pdbSecStr = removeBulges (pdbSecStr, 'H', 3)
        pdbSecStr = removeBulges (pdbSecStr, 'E', 3)

        if len(pdbSeq) == 0:
            raise Exception('empty pdb seq')
        if len(pdbSeq) != len(pdbSecStr):
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

        alignedpath += '_al'

        if not os.path.isfile (alignedpath):
            raise Exception ('alignment file %s not created, MSA error:\n%s' %
			     (alignedpath, p.stderr.read ()))

        aligned = parseFasta (open(alignedpath, 'r'))

        os.remove (alignedpath)
        os.remove (toalignpath)

        if aligned['template'].replace('-', '') != pdbSeq:
            raise Exception('MSA output mismatch:\npdbSeq:' + pdbSeq +
                            'aligned:' + aligned['template'])

        return aligned

aligner = alignService()
