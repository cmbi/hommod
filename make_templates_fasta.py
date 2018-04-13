import sys 
import os
import urllib2
import re
from argparse import ArgumentParser
import logging

settings = {}
filename = 'hommod/default_settings.py'
with open(filename) as config_file:
    exec(compile(config_file.read(), filename, 'exec'), settings)
env_settings = {}
filename = os.environ['HOMMOD_SETTINGS']
with open(filename) as config_file:
    exec(compile(config_file.read(), filename, 'exec'), env_settings)
settings.update(env_settings)
settings = {k:v for k, v in settings.iteritems() if k.isupper()}

from hommod.services.dssp import dssp
dssp.DSSP_DIR = settings['DSSP_DIR']

from hommod.controllers.blacklist import blacklister
blacklister.file_path = settings['BLACKLIST_FILE_PATH']

from hommod.controllers.fasta import write_fasta
from hommod.controllers.pdbfinder import parse_pdbfinder
from hommod.controllers.sequence import is_protein_sequence, is_nucleic_acid_sequence


_log = logging.getLogger(__name__)


P_UNALLOWED_HET_GROUPS = [
    re.compile('.*VANADATE.*'),
    re.compile('.*ARSENIC.*'),
    re.compile('.*URANIUM.*'),
    re.compile('NITRATE ION')
]

def contains_unallowed_het_groups(entry):
    het_groups = entry.find_one('HET-Groups')
    if het_groups is not None:
        for het_group in het_groups.children:
            name = het_group.find_one('Name')
            if name is not None:
                for p in P_UNALLOWED_HET_GROUPS:
                    if p.match(name.value):
                        _log.debug("{} contains unallowed {}".format(entry.find_one('ID').value, name.value))
                        return True
    return False


P_NUC = re.compile(r'[actgu\-]+')
P_PROT = re.compile(r'[ABCDEFGHIJKLMNOPQRSTUVWXYZ\-]+')

def get_sequences():
    sequences = {}

    for entry in parse_pdbfinder(settings['PDBFINDER2_FILE_PATH']):

        pdbid = entry.find_one('ID').value

        if blacklister.is_blacklisted(pdbid):
            _log.debug("{} is blacklisted".format(pdbid))
            continue

        if contains_unallowed_het_groups(entry):
            continue

        for chain in entry.find('Chain'):
            sequence = chain.find_one('Sequence')
            if sequence is not None:
                _log.debug("{} {} has sequence {}".format(pdbid, chain.value, sequence.value))

                if P_PROT.match(sequence.value):
                    fasta_key = 'pdb|%s|%s' % (pdbid.upper(), chain.value)

                    # PDBFINDER can contain multiple sequences with the same id.
                    # Always take the largest.
                    if fasta_key in sequences and len(sequences[fasta_key]) > len(sequence.value):
                        continue

                    # Blast cannot handle '-' in the sequence, so replace it with 'X'.
                    sequences[fasta_key] = sequence.value.replace('-', 'X')

                    _log.info("{} {}".format(fasta_key, sequences[fasta_key]))

    return sequences


if __name__ == "__main__":

    logging.basicConfig()
    if settings['DEBUG']:
        _log.setLevel(logging.DEBUG)

    parser = ArgumentParser(description='Make a fasta of all usable templates')
    parser.add_argument('output_file', help='the output fasta file')

    args = parser.parse_args()

    sequences = get_sequences()

    write_fasta(args.output_file, sequences)
