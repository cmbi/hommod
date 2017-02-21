from pathlib import Path

import logging

from celery import current_app as celery_app

from hommod_rest.services.modelutils import parseFasta, writeFasta
from hommod_rest.services.model import modeler
from hommod_rest.services.utils import (list_models_of, select_best_model,
                                        get_oldest_hg_sequence)

_log = logging.getLogger(__name__)


@celery_app.task()
def remodel_oldest_hg():
    fasta_path = get_oldest_hg_sequence()
    with open(fasta_path, 'r') as f:
        fasta = parseFasta(f)

    sequence = fasta.values()[0]

    # Touch the fasta first, setting the modification time to now.
    # This puts this fasta to the back of the sorted file list next time.
    Path(fasta_path).touch()

    modeler.modelProc(sequence, 'HUMAN', overwrite=True)


@celery_app.task()
def create_models_seq(sequence, species_id):
    _log.info("create_models_seq called with {} {}".format(sequence, species_id))

    modeler.modelProc(sequence, species_id)


@celery_app.task()
def create_model(sequence, species_id, residue_number, template_id):

    _log.info("create_model called with {} {} {} {}".format(sequence,
                                                            species_id,
                                                            residue_number,
                                                            template_id))
    if residue_number < 0 or residue_number > len(sequence):
        raise Exception("residue number %d out of sequence range" %
                        residue_number)

    paths = list_models_of(sequence, species_id, residue_number, template_id)

    _log.debug("create_model: {} present models for {} {} {} {}"
               .format(len(paths), sequence, species_id, residue_number,
                       template_id))

    if len(paths) <= 0:

        _log.debug("create_model: running modeler for {} {} {} {}"
                   .format(sequence, species_id, residue_number, template_id))

        paths = modeler.modelProc(sequence, species_id, residue_number, False,
                                  template_id)

    _log.debug("create_model: current models: %s" % str(paths))

    if len(paths) > 0:

        _log.debug("create_model: selecting best model for {} {} {} {} out of {}"
                   .format(sequence, species_id, residue_number, template_id, paths))

        best_model = select_best_model(sequence, species_id, residue_number,
                                       template_id)

        _log.debug("create_model: ending job with " +
                   "best model {} found for {} {} {} {}"
                   .format(best_model, sequence, species_id, residue_number,
                           template_id))
        return best_model
    else:
        # This isn't an error but I want to get an email when it happens so I
        # can keep an eye on it. There are too many warn log statements to set
        # the SMTP logger to send them.
        _log.error("create_model: 0 models created for  {} {} {} {}"
                   .format(sequence, species_id, residue_number, template_id))
        return ''
