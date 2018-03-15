#!/usr/bin/python

from flask import current_app as flask_app

from hommod_rest.services.modelutils import (parseFasta, getCoverageIdentity,
                                             idForSeq, TemplateID)
from hommod_rest.services.align import aligner
from hommod_rest.services.blast import blaster
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

# template: 3hrw
#     matched main target 7cd68612-3b7d-668d-756a-7f5573b87fa9 with chain B

ptargetchain = re.compile(
    r"^matched\s+(main\s+|)target\s+([A-Za-z0-9\-]+)\s+with chain\s+(\w)$")


def extract_info(tar_path):
    info = {'targets': {}}

    tf = tarfile.open(tar_path, 'r')

    infopath = os.path.join(
        os.path.splitext(
            os.path.basename(tar_path))[0], 'selected-targets.txt')

    _log.debug("extract info from %s" % infopath)

    if infopath in tf.getnames():
        for line in tf.extractfile(infopath):

            line = line.strip()
            if line.startswith('template:'):
                info['template'] = line.split(':')[1].strip()
            else:
                i_target = line.find(' target ')
                i_chain = line.find(' chain ')
                if i_target != -1 and i_chain != -1:

                    target = line[i_target:].split()[1]
                    chain = line[i_chain:].split()[1]

                    _log.debug("matched target chain line:\n\"%s\"\n" % line +
                               "chain=\'%s\', target=\'%s\'" % (chain, target))
                    info['targets'][chain] = target
    tf.close()

    _log.debug("got info {}".format(info))

    return info


def extract_template_id(tar_path):
    info = extract_info(tar_path)
    if 'template' in info and '_' in info['template']:
        try:
            pdbac, chain = info['template'].split('_')
        except:
            _log.debug("found no template in {}".format(tar_path))
            return None
        template_id = TemplateID(pdbac, chain)
        _log.debug("found template {} in {}".format(template_id, tar_path))
        return template_id
    else:
        _log.debug("found no template in {}".format(tar_path))
        return None


def extract_alignment(tar_path):

    _log.debug("extract alignment from %s" % tar_path)

    alignment = {}

    tf = tarfile.open(tar_path, 'r')

    fastapath = os.path.join(os.path.splitext(os.path.basename(tar_path))[0],
                             'align.fasta')
    if fastapath in tf.getnames():
        alignment = parseFasta(tf.extractfile(fastapath))

    tf.close()

    _log.debug("got model alignment {}".format(alignment))

    return alignment


def extract_model(tar_path):
    _log.debug("extract model from %s" % tar_path)

    contents = ''
    tf = tarfile.open(tar_path, 'r')

    pdbpath = os.path.join(os.path.splitext(os.path.basename(tar_path))[0],
                           'target.pdb')
    if pdbpath in tf.getnames():
        contents = tf.extractfile(pdbpath).read()

    tf.close()

    _log.debug("got model contents")

    return contents


def get_oldest_hg_sequence():
    fastadir = flask_app.config['HGFASTADIR']
    fastas = sorted([os.path.join(fastadir, filename)
        for filename in os.listdir(fastadir)], key=os.path.getmtime)
    return fastas[0]


def select_best_model(user_sequence, paths):
    bestID = 0.0
    best = None
    for path in paths:

        alignment = extract_alignment(path)

        templateID = None
        for ID in alignment.keys():
            if ID != 'target':
                templateID = ID
                break

        if templateID is None:
            raise Exception("no template found in alignment of %s" % path)

        main_i = -1
        targetseqs = alignment['target'].split('|')
        templateseqs = alignment[templateID].split('|')
        for i in range(len(targetseqs)):

            # target sequence might have been slightly altered during the procedure
            # That's why we must align here.
            aligned = aligner.clustal_align({
                't': targetseqs[i].replace('-', ''),
                'm': user_sequence})
            pcov, pid = getCoverageIdentity(aligned['t'], aligned['m'])
            if pid > 95.0:

                main_i = i
                break

        if main_i == -1:
            raise Exception("main target sequence not found in " + path + " alignment:\n" +
                            user_sequence + " should match one of " +
                            str(targetseqs))

        pcov, pid = \
            getCoverageIdentity(targetseqs[main_i], templateseqs[main_i])

        if pid > bestID:
            best = path
            bestID = pid

    return best


def blast_models(sequence, species, position, template_id):

    approved = []

    hits = blaster.blast_models(sequence)
    for hitID in hits:

        path, chain = hitID.split('|')

        name = os.path.splitext(os.path.basename(path))[0]

        if ('_' + species + '_') not in name:
            continue

        if not os.path.isfile(path):
            continue

        info = extract_info(path)

        if template_id is not None and info['template'].split('_')[0] != template_id.pdbac:
            continue

        for alignment in hits[hitID]:
            if position < alignment.querystart or position > alignment.queryend:
                continue

            length = alignment.queryend - alignment.querystart
            if length < 20:
                continue

            # must be 100% identity
            if alignment.queryalignment != alignment.subjectalignment:
                continue

            approved.append(path)
            break

    return approved


def get_model_path(model_id):
    path = os.path.join(flask_app.config['MODELDIR'], '%s.tgz' % model_id)
    hg_path = os.path.join(flask_app.config['HGMODELDIR'], '%s.tgz' % model_id)

    _log.debug("checking path {}".format(path))
    _log.debug("checking path {}".format(hg_path))

    if os.path.isfile(hg_path):

        _log.debug("{} exists".format(hg_path))

        if os.path.isfile(path):
            _log.debug("{} exists".format(path))

            # Return latest file path
            if os.path.getmtime(path) < os.path.getmtime(hg_path):
                return hg_path
            else:
                return path
        else:
            return hg_path
    elif os.path.isfile(path):
        _log.debug("{} exists".format(path))

        return path
    else:
        return None


def list_models_of(sequence, species, position, template_id):

    _log.debug("list_models_of {} {}".format(sequence, species))

    seq_id = idForSeq(sequence)

    _log.debug("list_models_of {} {}".format(seq_id, species))

    wildcard = os.path.join(flask_app.config['MODELDIR'],
                            '%s_%s_*-*.tgz' % (seq_id, species))
    hg_wildcard = os.path.join(flask_app.config['HGMODELDIR'],
                               '%s_%s_*-*.tgz' % (seq_id, species))

    _log.debug("looking for " + wildcard + " and " + hg_wildcard)

    l = []
    for f in (glob(wildcard) + glob(hg_wildcard)):

        _log.debug("inspecting {}".format(f))

        age = time() - os.path.getmtime(f)
        if age >= flask_app.config['MAX_MODEL_DAYS'] * 24 * 60 * 60:
            continue

        name = os.path.splitext(os.path.basename(f))[0]

        _log.debug("name is {}".format(name))

        if template_id is not None:
            _log.debug("checking whether {} is the template of {}".format(
                template_id, f))
            f_template_id = extract_template_id(f)
            if f_template_id != template_id:
                _log.debug("skipping {}, because it does not have {} as template, but {}"
                           .format(f, template_id, f_template_id))
                continue

        if position is None:
            _log.debug("add {} to list".format(f))
            l.append(f)
            continue

        name_params = name.split('_')
        ran = name_params[2].split('-')

        start = int(ran[0])
        end = int(ran[1])

        if position >= start and position <= end:
            try:
                model = extract_model(f)
                alignment = extract_alignment(f)
            except:
                _log.error("couldn't get model or alignment from %s" % f)
                continue

            if len(model) == 0:
                _log.error("empty model in %s" % f)
                continue

            inputtargetseq = sequence[start - 1: end]
            modeltargetseqs = \
                alignment['target'].replace('-', '').replace('.', '').split('|')

            fragment_present = False
            for modeltargetseq in modeltargetseqs:
                if inputtargetseq in modeltargetseq:
                    fragment_present = True
                    break

            if not fragment_present:
                _log.error(
                    "sequence mismatch in %s:\nquery target: %s\naligned target: %s"
                    % (f, inputtargetseq, alignment['target']))
                continue

            _log.debug("add {} to list".format(f))
            l.append(f)

    _log.debug("found models {}".format(l))
    return l
