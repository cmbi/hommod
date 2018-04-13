import os
import logging
import traceback

from filelock import FileLock
from celery import current_app as celery_app
from celery import group
from celery.signals import task_failure

from hommod.controllers.model import modeler
from hommod.controllers.storage import model_storage
from hommod.controllers.domain import domain_aligner
from hommod.models.error import InitError


_log = logging.getLogger(__name__)


@celery_app.task()
def create_models(target_sequence, target_species_id, require_resnum=None, chosen_template_id=None):

    target_species_id = target_species_id.upper()

    sequence_id = model_storage.get_sequence_id(target_sequence)
    lock_name = "lock_search_%s_%s_%s_%s" % (sequence_id,
                                             target_species_id,
                                             str(require_resnum),
                                             str(chosen_template_id))

    if model_storage.model_dir is None:
        raise InitError("model directory is not set")

    lock_path = os.path.join(model_storage.model_dir, lock_name)
    with FileLock(lock_path):

        model_paths = model_storage.list_models(sequence_id, target_species_id,
                                                require_resnum, chosen_template_id)
        if len(model_paths) > 0:
            return model_paths
        else:
            overlapping_domain_alignments = \
                domain_aligner.get_domain_alignments(target_sequence,
                                                     require_resnum,
                                                     chosen_template_id)

            job = group([create_model.subtask((target_sequence,
                                               target_species_id,
                                               domain_alignment))
                       for domain_alignment in overlapping_domain_alignments])

            result = job.apply_async()

            return result.join()


@celery_app.task()
def create_model(target_sequence, target_species_id, domain_alignment):

    target_species_id = target_species_id.upper()

    return modeler.build_model(target_sequence, target_species_id, domain_alignment)


@task_failure.connect
def task_failure_handler(task_id, exception, *args, **kwargs):
    _log.error(traceback.format_tb(kwargs['traceback']))
