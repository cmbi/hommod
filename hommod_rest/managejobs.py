from hashlib import md5
from flask import current_app as flask_app


def add_create_model_job(sequence, species_id, position):

    jobid = 'create_model_%s_%s_%i' % (md5(sequence).hexdigest(), species_id,
                                       position)

    # See if the job with this ID is already running/completed:
    from hommod_rest.tasks import create_model
    prevResult = create_model.AsyncResult(jobid)
    if prevResult.state != 'PENDING' and \
        not (prevResult.state == 'FAILURE' and
             flask_app.config['RETRY_FAILURE']):

        return prevResult.task_id

    result = create_model.apply_async((sequence, species_id, position),
                                      task_id=jobid)

    # Need to set it to something other than PENDING !
    create_model.update_state(result.task_id, 'SENT')

    return result.task_id
