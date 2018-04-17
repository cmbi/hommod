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

from hommod.controllers.pdb import parse_seqres_from_string
from hommod.controllers.fasta import write_fasta

from hommod.controllers.storage import model_storage
model_storage.model_dir = settings['MODEL_DIR']


_log = logging.getLogger(__name__)


def get_sequences():

    sequences = {}

    for tar_path in model_storage.list_all_models():
        try:
            contents = model_storage.extract_model(tar_path)
        except:
            continue

        seqres = parse_seqres_from_string(contents)
        for chain_id in seqres:
            sequences[tar_path + '|' + chain_id] = ''.join([aa.letter for aa in seqres[chain_id]])

    return sequences


if __name__ == "__main__":

    logging.basicConfig()
    if settings['DEBUG']:
        _log.setLevel(logging.DEBUG)

    parser = ArgumentParser(description='Make a fasta of all usable models')
    parser.add_argument('output_file', help='the output fasta file')

    args = parser.parse_args()

    sequences = get_sequences()

    write_fasta(args.output_file, sequences)
