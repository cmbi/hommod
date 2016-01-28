import logging
import inspect
import re

from flask import Blueprint, jsonify, request, Response, render_template

from hommod_rest.services.utils import (extract_info, extract_alignment,
                                        extract_model, list_models_of)

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

    _log.info("update cache request for ( sequence: %s, species: %s )"
               % (sequence, species_id))

    if not (sequence and species_id):
        return jsonify({'error': 'invalid input'}), 400

    from hommod_rest.tasks import create_models_seq
    result = create_models_seq.apply_async((sequence, species_id))

    _log.info("cache will update in job %s" % result.task_id)

    return jsonify({'jobid': result.task_id})

@bp.route('/has_model/', methods=['POST'])
def has_model ():

    """
    Ask whether the given model already exists or not.

    :param sequence: target sequence for the model
    :param position: position of the required residue in the sequence, starting 1
    :param species_id: uniprot species id for the model
    :return: True if the model is already there, false otherwise
    """

    sequence = request.form.get('sequence', None)
    position = request.form.get('position', None)
    species_id = request.form.get('species_id', None)

    _log.info ("has_model request recieved for ( sequence: %s, species: %s, position: %s )"
               % (sequence, species_id, position))

    paths = list_models_of (sequence, species_id, position)

    if len (paths) > 0:

        _log.info ("has_model: returning true for %s" % str (paths))
        return True
    else:
        _log.info ("has_model: returning false for %s" % str (paths))
        return False

@bp.route('/submit/', methods=['POST'])
def submit ():

    """
    Request a model for the given target 'sequence', residue position
    and 'species_id'.

    :param sequence: target sequence for the model
    :param position: position of the required residue in the sequence, starting 1
    :param species_id: uniprot species id for the model
    :return: a json object, containing the field 'jobid'
    """

    sequence = request.form.get('sequence', None)
    position = request.form.get('position', None)
    species_id = request.form.get('species_id', None)

    _log.info ("submit request recieved for ( sequence: %s, species: %s, position: %s )"
               % (sequence, species_id, position))

    if not (sequence and position and species_id):

        _log.warn ("submit request did not contain all required input data");

        return jsonify({'error': 'invalid input'}), 400
    try:
        position = int(position)
    except:
        _log.warn ("submit request did not contain an integer position");

        return jsonify({'error': 'expected integer for position'}), 400

    _log.debug("submitted ( sequence: %s, species: %s, position: %i )"
               % (sequence, species_id, position))

    from hommod_rest.tasks import create_model
    result = create_model.apply_async((sequence, species_id, position))

    _log.info ("created job %s" % result.task_id)

    return jsonify({'jobid': result.task_id})


@bp.route('/status/<jobid>/', methods=['GET'])
def status (jobid):

    """
    Request the status of a job.

    :param jobid: the jobid returned by 'submit'
    :return: Either PENDING, STARTED, SUCCESS, FAILURE, RETRY, or REVOKED.
    """

    _log.info ("status request for job %s" % jobid)

    from hommod_rest.application import celery
    result = celery.AsyncResult(jobid)

    response = {'status': result.status}
    if result.failed ():
        response ['message'] = result.traceback

    _log.info ("found status %s for job %s" % (result.status, jobid))

    return jsonify (response)


@bp.route('/get_model_file/<jobid>.pdb', methods=['GET'])
@bp.route('/get_model_file/<jobid>.PDB', methods=['GET'])
def get_model_file (jobid):

    """
    Get the pdb file, created by the modeling job.

    :param jobid: the jobid returned by 'submit'
    :return: The pdb file created by the job. If the job status is not SUCCESS, this method returns an error.
    """

    _log.info ("model request for job %s" % jobid)

    from hommod_rest.application import celery
    result = celery.AsyncResult(jobid)
    path = result.result
    if not path:  # no model could be created

        _log.warn('no model was created for job %s' % jobid)

#        return '', 404
        return ''

    try:
        contents = extract_model(path)
    except:
        _log.warn('failed to get all data from %s' % path)
        return '', 500

    _log.info ("model successfully retrieved for job %s" % jobid)

    return Response(contents, mimetype='chemical/x-pdb')


@bp.route('/get_metadata/<jobid>/', methods=['GET'])
def get_metadata (jobid):

    """
    Get the metadata of the model, created by the modeling job.

    :param jobid: the jobid returned by 'submit'
    :return: The json metadata object. If the job status is not SUCCESS, this method returns an error.
    """

    _log.info ("metadata request for job %s" % jobid)

    from hommod_rest.application import celery
    result = celery.AsyncResult(jobid)
    path = result.result
    if not path:
#        return {'error': 'no model could be created'}, 404
	    return {}

    try:
        data = extract_info(path)
        data['alignment'] = extract_alignment(path)
    except:
        _log.warn ('failed to get all data from %s' % path)
        return 'data not available', 500

    _log.info ("metadata successfully retrieved for job %s" % jobid)

    return jsonify(data)


@bp.route ('/')
def docs ():

    _log.info ("endpoint request for hommod docs")

    fs = [update_cache, submit, has_model, status, get_model_file, get_metadata]
    docs = {}
    for f in fs:
        src = inspect.getsourcelines (f)
        m = re.search(r"@bp\.route\('([\w\/\<\>\.]*)', methods=\['([A-Z]*)']\)",
                      src[0][0])

        if not m:  # pragma: no cover

            _log.warn ("Unable to document function '{}'".format(f))
            continue

        url = m.group(1)
        method = m.group(2)
        docstring = inspect.getdoc(f)
        docs[url] = (method, docstring)

    _log.info ("rendering template for hommod docs")

    return render_template('api/docs.html', docs=docs)

