import logging
import inspect
import re

from flask import Blueprint, jsonify, request, Response

from hommod_rest.services.utils import (extract_info, extract_alignment,
                                        extract_model)

_log = logging.getLogger(__name__)


bp = Blueprint('hommod', __name__, url_prefix='/api')


@bp.route('/update_cache/', methods=['POST'])
def update_cache ():

    """
    Request a rebuild for all models for the given target 'sequence'
    and 'species_id' that could possibly be built.

    :param sequence: target sequence for the model
    :param species_id: uniprot species id for the model
    :return: a json object, containing the field 'jobid'
    """

    sequence = request.form.get('sequence', None)
    species_id = request.form.get('species_id', None)

    if not (sequence and species_id):
        return jsonify({'error': 'invalid input'}), 400

    _log.debug("update resquest for ( sequence: %s, species: %s )"
               % (sequence, species_id))

    from hommod_rest.tasks import create_models_seq
    result = create_models_seq.apply_async((sequence, species_id))

    return jsonify({'jobid': result.task_id})

@bp.route('/submit/', methods=['POST'])
def submit ():
    sequence = request.form.get('sequence', None)
    position = request.form.get('position', None)
    species_id = request.form.get('species_id', None)

    if not (sequence and position and species_id):
        return jsonify({'error': 'invalid input'}), 400
    try:
        position = int(position)
    except:
        return jsonify({'error': 'expected integer for position'}), 400

    _log.debug("submitted ( sequence: %s, species: %s, position: %i )"
               % (sequence, species_id, position))

    from hommod_rest.tasks import create_model
    result = create_model.apply_async((sequence, species_id, position))

    return jsonify({'jobid': result.task_id})


@bp.route('/status/<jobid>/', methods=['GET'])
def status (jobid):

    """
    Request the status of a job.

    :param jobid: the jobid returned by 'submit'
    :return: Either PENDING, STARTED, SUCCESS, FAILURE, RETRY, or REVOKED.
    """

    from hommod_rest.application import celery
    result = celery.AsyncResult(jobid)
    return jsonify({'status': result.status})


@bp.route('/get_model_file/<jobid>.pdb', methods=['GET'])
def get_model_file (jobid):

    """
    Get the pdb file, created by the modeling job.

    :param jobid: the jobid returned by 'submit'
    :return: The pdb file created by the job. If the job status is not SUCCESS, this method returns an error.
    """

    from hommod_rest.application import celery
    result = celery.AsyncResult(jobid)
    path = result.result
    if not path:  # no model could be created
        return '', 404

    try:
        contents = extract_model(path)
    except:
        _log.warn('failed to get all data from %s' % path)
        return '', 500

    return Response(contents, mimetype='chemical/x-pdb')


@bp.route('/get_metadata/<jobid>/', methods=['GET'])
def get_metadata (jobid):

    """
    Get the metadata of the model, created by the modeling job.

    :param jobid: the jobid returned by 'submit'
    :return: The json metadata object. If the job status is not SUCCESS, this method returns an error.
    """

    from hommod_rest.application import celery
    result = celery.AsyncResult(jobid)
    path = result.result
    if not path:
        return {'error': 'no model could be created'}, 404

    try:
        data = extract_info(path)
        data['alignment'] = extract_alignment(path)
    except:
        _log.warn('failed to get all data from %s' % path)
        return 'data not available', 500

    return jsonify(data)


@bp.route ('/')
def docs ():

    p = re.compile (r"\@bp\.route\s*\(\'([\w\/\<\>]*)\'\)")

    fs = [annotations, entries]
    docs = {}
    for f in fs:
        src = inspect.getsourcelines (f)
        m = p.search (src[0][0])
        if not m:  # pragma: no cover
            _log.debug("Unable to document function '{}'".format(f))
            continue

        url = m.group(1)
                                                                                                                                                              55,1          Top
