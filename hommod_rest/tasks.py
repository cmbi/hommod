import logging

from celery import current_app as celery_app

from hommod_rest.services.model import modeler
from hommod_rest.services.utils import list_models_of, select_best_model

_log = logging.getLogger(__name__)



@celery_app.task()
def create_models_seq(sequence, species_id):

    _log.info("creating models for {} {}".format(sequence, species_id))

    modeler.modelProc(sequence, species_id)

@celery_app.task()
def create_model(sequence, species_id, residue_number):

    _log.info("creating model with {} {} {}".format(sequence, species_id,
                                                    residue_number))

    paths = list_models_of(sequence, species_id, residue_number)

    _log.debug("{} models already there for {} {} {}"
               .format(len(paths), sequence, species_id, residue_number))

    if len(paths) <= 0:
        paths = modeler.modelProc(sequence, species_id, residue_number, False)

    if len(paths) > 0:
        return select_best_model(sequence, species_id, residue_number)
    else:
        _log.warn("failed a model: {} {} {}"
                  .format(sequence, species_id, residue_number))
        raise Exception("unable to create model for {} {} {}"
                        .format(sequence, species_id, residue_number))
