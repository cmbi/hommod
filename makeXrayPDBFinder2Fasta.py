#!/usr/bin/python

import sys
import os
import urllib2
import re

settings = {}
filename = 'hommod_rest/default_settings.py'
with open(filename) as config_file:
    exec(compile(config_file.read(), filename, 'exec'), settings)
env_settings = {}
filename = os.environ['HOMMOD_REST_SETTINGS']
with open(filename) as config_file:
    exec(compile(config_file.read(), filename, 'exec'), env_settings)
settings.update(env_settings)
settings = {k:v for k, v in settings.iteritems() if k.isupper()}


def has_dssp(pdbid):
    return os.path.isfile(os.path.join(settings['DSSPDIR'], '%s.dssp' % pdbid))


def writeFasta(seqs, path):
    f = open(path, 'w')
    for key in seqs:
        f.write('>%s\n%s\n' % (key, seqs[key]))
    f.close()


def gather(blacklisted_templates):
    # First gather all sequences from the pdbFinder:
    pNuc = re.compile(r'[actgu\-]+')
    pProt = re.compile(r'[ABCDEFGHIJKLMNOPQRSTUVWXYZ\-]+')

    seqs = {}
    chain_id = None
    pdbid = None
    section = None
    method = None
    resolution = None

    # Parse entire set of available sequences:
    with open(settings['PDBFINDER2'], 'r') as source:
        for line in source:
            if line.startswith('//'):  # comment
                continue
            i = line.find(':')
            if i < 0:  # no colon
                continue
            line_key = line[:i].strip()
            line_value = line[i + 1:].strip()

            if not line[0].isspace():
                section = line_key.strip()  # Chain, HET-Groups, Exp-Method, etc.

            if line_key == 'ID':
                pdbid = line_value.lower()

            elif line_key == 'Chain':
                chain_id = line_value

            elif line_key == 'Exp-Method':
                method = line_value

            elif line_key == 'Resolution' and section == 'Exp-Method':
                resolution = float(line_value)

            elif line_key == 'Sequence':
                seq = line_value

                if len(seq) > 0 and method == 'X' and has_dssp(pdbid) and \
                        not pdbid in blacklisted_templates and \
                        pProt.match(seq) and not pNuc.match(seq):

                    # Generate key, unique for pdb file and chain id:
                    fasta_key = 'pdb|%s|%s' % (pdbid.upper(), chain_id)

                    # Take the longest sequence in the chain:
                    if fasta_key not in seqs or \
                            len(seqs[fasta_key]) < len(seq):
                        seqs[fasta_key] = seq.replace('-', 'X')
        source.close()

    return seqs


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print 'Usage: %s [output file]' % sys.argv[0]
        sys.exit(0)

    out_file = sys.argv[1]

    blacklisted_templates = []
    with open(settings['TEMPLATE_BLACKLIST'], 'r') as f:
        blacklisted_templates = f.read().lower().split()

    seqs = gather(blacklisted_templates)
    writeFasta(seqs, out_file)
