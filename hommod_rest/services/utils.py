#!/usr/bin/python

from flask import current_app as flask_app

from hommod_rest.services.modelutils import (parseFasta, getCoverageIdentity,
                                             idForSeq)
from hommod_rest.services.align import aligner
from time import time

import re
import os
import tarfile
from glob import glob

import logging

_log = logging.getLogger(__name__)
sh = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
sh.setFormatter(formatter)
_log.addHandler(sh)
_log.setLevel(logging.DEBUG)

ptargetchain = re.compile(
    r"^modeling\s+(main\s+|)target\s+([A-Za-z0-9\-]+)\s+on chain\s+(\w)$")


def extract_info(tar_path):

    info = {'targets': {}}

    tf = tarfile.open(tar_path, 'r')

    infopath = os.path.join(
        os.path.splitext(os.path.basename(tar_path))[0], 'selected-targets.txt')
    if infopath in tf.getnames():
        for line in tf.extractfile(infopath):

            line = line.strip()
            if line.startswith('template:'):
                info['template'] = line.split(':')[1].strip()
            else:
                m = ptargetchain.match(line)
                if m:
                    chain = m.group(3)
                    target = m.group(2)
                    info['targets'][chain] = target

    tf.close()

    return info


def extract_alignment(tar_path):

    alignment = {}

    tf = tarfile.open(tar_path, 'r')

    fastapath = os.path.join(os.path.splitext(os.path.basename(tar_path))[0],
                             'align.fasta')
    if fastapath in tf.getnames():
        alignment = parseFasta(tf.extractfile(fastapath))

    tf.close()

    return alignment


def extract_model(tar_path):

    contents = ''
    tf = tarfile.open(tar_path, 'r')

    pdbpath = os.path.join(os.path.splitext(os.path.basename(tar_path))[0],
                           'target.pdb')
    if pdbpath in tf.getnames():
        contents = tf.extractfile(pdbpath).read()

    tf.close()

    return contents


def select_best_model(sequence, species, position):

    bestID = 0.0
    best = None
    for path in list_models_of (sequence, species, position):

        alignment = extract_alignment(path)

        templateID = None
        for ID in alignment.keys():
            if ID != 'target':
                templateID = ID
                break

        if templateID is None:
            raise Exception ("no template found in alignment of %s" % path)

        main_i = -1
        targetseqs = alignment['target'].split('|')
        templateseqs = alignment[templateID].split('|')
        for i in range(len(targetseqs)):

            # target sequence might have been slightly altered during the procedure
            # That's why we must align here.
            aligned = aligner.clustalAlign ({'t':targetseqs[i].replace('-', ''),'m':sequence})
            pcov, pid = getCoverageIdentity (aligned ['t'], aligned ['m'])
            if pid > 95.0:

                main_i = i
                break

        if main_i == -1:
            raise Exception ("main target sequence not found in alignment:\n" +
                            sequence + " should match one of " +
                            str(targetseqs))

        pcov, pid = \
            getCoverageIdentity(targetseqs[main_i], templateseqs[main_i])

        if pid > bestID:
            best = path
            bestID = pid

    _log.debug ("best model for %s %s %i is %s" % (idForSeq(sequence), species, position, best))
    return best


def list_models_of(sequence, species, position):

    h = idForSeq(sequence)

    wildcard = os.path.join(flask_app.config['MODELDIR'],
                            '%s_%s_*-*.tgz' % (h, species))
    _log.debug("looking for " + wildcard)

    l = []
    for f in glob(wildcard):

        age = time() - os.path.getmtime(f)
        if age >= flask_app.config['MAX_MODEL_DAYS']*24*60*60:
            continue

        name = os.path.splitext(os.path.basename(f))[0]
        ran = name.split('_')[-1].split('-')

        start = int(ran[0])
        end = int(ran[1])

        if position >= start and position <= end:

            try:
                model = extract_model (f)
                alignment = extract_alignment (f)
            except:
                _log.error ("couldn't get model or alignment from %s" % f)
                continue

            if len (model) == 0:
                _log.error ("empty model in %s" % f)
                continue

            inputtargetseq = sequence [start - 1: end]
            modeltargetseq = alignment ['target'].replace ('-', '').replace ('.', '')
            if modeltargetseq != inputtargetseq:
                _log.error ("sequence mismatch in %s:\ninput: %s\nalign: %s" % (f, inputtargetseq, modeltargetseq))

            _log.debug ("add %s to list of models for %s %s %i" % (f, h, species, position))
            l.append(f)

    _log.debug ("found %i models for %s %s %i" % (len (l), h, species, position))
    return l
