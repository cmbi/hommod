import logging

from celery import current_app as celery_app

from hommod_rest.services.model import modeler
from hommod_rest.services.utils import list_models_of

_log = logging.getLogger(__name__)


@celery_app.task()
def create_model(sequence, species_id, residue_number):

    _log.info("creating model with {} {} {}".format(sequence, species_id,
                                                    residue_number))

    paths = []
    paths = list_models_of(sequence, species_id, residue_number)

    if len(paths) <= 0:
        paths = modeler.modelProc(sequence, species_id, residue_number, False)

    return paths
