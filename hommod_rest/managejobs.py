# from hashlib import md5
# from flask import current_app as flask_app


def create_model(sequence, species_id, position):

    # jobid = 'create_model_%s_%s_%i' % (md5(sequence).hexdigest(), species_id,
    #                                    position)

    # See if the job with this ID is already running/completed:
    # from hommod_rest.application import celery
    # prevResult = celery.AsyncResult(jobid)
    # if prevResult.state != 'PENDING' and \
    #     not (prevResult.state == 'FAILURE' and
    #          flask_app.config['RETRY_FAILURE']):
    #     return prevResult.task_id

    from hommod_rest.tasks import create_model
    result = create_model.apply_async((sequence, species_id, position))
    # result = create_model.apply_async((sequence, species_id, position),
    #                                   task_id=jobid)

    return result.task_id
