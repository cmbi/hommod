#!/usr/bin/python

from flask import current_app as flask_app

from hommod_rest.services.modelutils import (parseFasta, getCoverageIdentity,
                                             idForSeq, TemplateID)
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

#template: 3hrw
#    matched main target 7cd68612-3b7d-668d-756a-7f5573b87fa9 with chain B

ptargetchain = re.compile(
    r"^matched\s+(main\s+|)target\s+([A-Za-z0-9\-]+)\s+with chain\s+(\w)$")


def extract_info(tar_path):

    info = {'targets': {}}

    tf = tarfile.open(tar_path, 'r')


    infopath = os.path.join(
                 os.path.splitext(
                   os.path.basename(tar_path))[0], 'selected-targets.txt')

    _log.debug ("extract info from %s" % infopath)

    if infopath in tf.getnames():
        for line in tf.extractfile(infopath):

            line = line.strip()
            if line.startswith('template:'):
                info['template'] = line.split(':')[1].strip()
            else:
                i_target = line.find (' target ')
                i_chain = line.find (' chain ')
                if i_target != -1 and i_chain != -1:

                    target = line [i_target:].split () [1]
                    chain = line [i_chain:].split () [1]

                    _log.debug ("matched target chain line:\n\"%s\"\n" % line +
                                "chain=\'%s\', target=\'%s\'" % (chain, target))
                    info['targets'][chain] = target

    tf.close()

    return info


def extract_template_id(tar_path):

    info = extract_info(tar_path)
    if 'template' in info and '_' in info['template']:
        pdbac, chain = info['template'].split('_')
        template_id = TemplateID(pdbac, chain)
        _log.debug("found template {} in {}".format(template_id, tar_path))
        return template_id
    else:
        _log.debug("found no template in {}".format(tar_path))
        return None


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


def select_best_model(sequence, species, position, template):

    bestID = 0.0
    best = None
    for path in list_models_of (sequence, species, position, template):

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

    _log.debug ("best model for %s %s %i %s is %s" % (idForSeq(sequence),
                species, position, str(template), best))
    return best


def list_models_of(sequence, species, position, template_id):

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

        if template_id:
            _log.debug("checking whether {} is the template of {}".format(
                        template_id, f))
            f_template_id = extract_template_id(f)
            if f_template_id != template_id:
                _log.debug("skipping {}, because it does not have {} as template, but {}"
                            .format(f, template_id, f_template_id))
                continue

        name_params = name.split('_')
        ran = name_params[2].split('-')

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
            modeltargetseqs = alignment ['target'].replace ('-', '').replace ('.', '').split ('|')

            fragment_present = False
            for modeltargetseq in modeltargetseqs:
                if inputtargetseq in modeltargetseq:
                    fragment_present = True
                    break

            if not fragment_present:
                _log.error ("sequence mismatch in %s:\nquery target: %s\naligned target: %s" % (f, inputtargetseq, alignment ['target']))
                continue

            _log.debug ("add %s to list of models for %s %s %i" % (f, h, species, position))
            l.append (f)

    _log.debug ("found %i models for %s %s %i" % (len (l), h, species, position))
    return l
