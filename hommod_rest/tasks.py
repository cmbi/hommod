import logging

from celery import current_app as celery_app

from hommod_rest.services.model import modeler
from hommod_rest.services.utils import list_models_of, select_best_model

_log = logging.getLogger(__name__)

@celery_app.task()
def create_models_seq(sequence, species_id):

    _log.info("create_models_seq called with {} {}".format(sequence, species_id))

    modeler.modelProc(sequence, species_id)

@celery_app.task()
def create_model(sequence, species_id, residue_number):

    _log.info("create_model called with {} {} {}".format(sequence, species_id,
                                                    residue_number))

    paths = list_models_of(sequence, species_id, residue_number)

    _log.debug("create_model: {} present models for {} {} {}"
               .format(len(paths), sequence, species_id, residue_number))

    if len(paths) <= 0:

        _log.debug ("create_model: running modeler for {} {} {}"
                    .format(sequence, species_id, residue_number))

        paths = modeler.modelProc (sequence, species_id, residue_number, False)

    _log.debug("create_model: current models: %s" % str (paths))

    if len(paths) > 0:

        _log.debug("create_model: selecting best model for {} {} {} out of {}"
                   .format (sequence, species_id, residue_number, paths))

        best_model = select_best_model (sequence, species_id, residue_number)

        _log.debug ("create_model: ending job with " + 
                    "best model {} found for {} {} {}"
                    .format (best_model, sequence, species_id, residue_number))
        return best_model
    else:
        _log.warn("create_model: 0 models created for  {} {} {}"
                  .format(sequence, species_id, residue_number))
	return ''
