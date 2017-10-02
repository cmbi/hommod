#!/usr/bin/python

import sys
import os
import tarfile
from glob import glob


def parse_fasta(file_):
    seqs = {}
    currentID = None
    for line in file_:
        if line.startswith('>'):
            currentID = line[1:].strip()
            seqs[currentID] = ''
        else:
            seqs[currentID] += line.strip()
    return seqs


dict_aa3to1 = {
    'ALA': 'A',
    'CYS': 'C',
    'ASP': 'D',
    'GLU': 'E',
    'PHE': 'F',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LYS': 'K',
    'LEU': 'L',
    'MET': 'M',
    'ASN': 'N',
    'PYL': 'O',
    'PRO': 'P',
    'ARG': 'R',
    'SER': 'S',
    'THR': 'T',
    'SEC': 'U',
    'VAL': 'V',
    'TRP': 'W',
    'TYR': 'Y'
}

def aa3to1(code):
    if code in dict_aa3to1:
        return dict_aa3to1[code]
    else:
        return 'X'


def get_sequences(tar_path):
    name = os.path.splitext(os.path.basename(tar_path))[0]

    seqs = {}
    with tarfile.open(tar_path, 'r') as tf:
        if '%s/target.pdb' % name not in tf.getnames():
            return {}

        f = tf.extractfile('%s/target.pdb' % name)
        for line in f:
            line = line.decode("utf-8")

            if line.startswith('SEQRES'):
                chain = line[11]
                if chain not in seqs:
                    seqs[chain] = ''

                for code in line[19:].split():
                    seqs[chain] += aa3to1(code)
        f.close()

    return seqs


if len(sys.argv) != 3:
    print("Usage: %s [models dir] [output fasta]" % sys.argv[0])
    sys.exit(1)

with open(sys.argv[2], 'w') as f:
    for path in glob('%s/*.tgz' % sys.argv[1]):
        name = os.path.splitext(os.path.basename(path))[0]

        seqs = get_sequences(path)
        for chain in seqs:
            f.write('>%s_%s\n%s\n' % (name, chain, seqs[chain]))

