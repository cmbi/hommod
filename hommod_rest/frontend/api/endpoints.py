import logging

from flask import Blueprint, jsonify, request, Response, render_template

from hommod_rest.services.utils import(extract_info, extract_alignment,
                                       extract_model, list_models_of)
from hommod_rest.services.modelutils import TemplateID

_log = logging.getLogger(__name__)


bp = Blueprint('hommod', __name__, url_prefix='/api')

@bp.route('/update_cache/', methods=['POST'])
def update_cache():

    """
    Request a rebuild for all models for the given target 'sequence'
    and 'species_id' that could possibly be built.

    :param sequence: target sequence for the model
    :param species_id: uniprot species id for the model
    :return: a json object, containing the field 'jobid'
    """

    sequence = request.form.get('sequence', None)
    species_id = request.form.get('species_id', None)

    _log.info("endpoints.update_cache request for( sequence: %s, species: %s )"
               %(sequence, species_id))

    if not(sequence and species_id):
        return jsonify({'error': 'invalid input'}), 400

    from hommod_rest.tasks import create_models_seq
    result = create_models_seq.apply_async((sequence, species_id))

    _log.debug("endpoints: cache will update in job %s" % result.task_id)

    return jsonify({'jobid': result.task_id})

@bp.route('/has_model/', methods=['POST'])
def has_model():

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
    template_id = request.form.get('template_id', None)

    try:
        position = int(position)
    except:
        _log.error("endpoints.submit: submit request did not contain an integer position");

        return jsonify({'error': 'expected integer for position'}), 400

    if type(template_id) == str and '-' in template_id:
        ac, chain = template_id.split('-')
        template_id = TemplateID(ac, chain)
    else:
        template_id = None

    _log.info("endpoints.has_model request recieved for( " +
              "sequence: %s, species: %s, position: %s ,template %s)"
              %(sequence, species_id, position, template_id))

    paths = list_models_of(sequence, species_id, position, template_id)

    if len(paths) > 0:

        _log.debug("endpoints.has_model: returning true for %s" % str(paths))
        return True
    else:
        _log.debug("endpoints.has_model: returning false for %s" % str(paths))
        return False

@bp.route('/submit/', methods=['POST'])
def submit():

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
    template_id = request.form.get('template_id', None)

    _log.info(("endpoints.submit: request recieved for( " +
               "sequence: %s, species: %s, position: %s, template: %s)")
               %(sequence, species_id, position, template_id))

    if not(sequence and position and species_id):

        _log.error("endpoints.submit: submit request did not contain all required input data");

        return jsonify({'error': 'invalid input'}), 400

    species_id = species_id.upper()
    try:
        position = int(position)
    except:
        _log.error("endpoints.submit: submit request did not contain an integer position");

        return jsonify({'error': 'expected integer for position'}), 400

    if '_' in template_id:
        ac, chain = template_id.split('_')
        template_id = TemplateID(ac, chain)
    else:
        template_id = None

    _log.debug("endpoints.submit: submitted ( sequence: %s, species: %s, position: %i, template: %s)"
               %(sequence, species_id, position, str(template_id)))

    from hommod_rest.tasks import create_model
    result = create_model.apply_async((sequence, species_id, position,
                                       template_id))

    _log.debug("endpoints.submit: created job %s, current_status=%s" % 
               (result.task_id, result.status))

    return jsonify({'jobid': result.task_id})


@bp.route('/status/<jobid>/', methods=['GET'])
def status(jobid):

    """
    Request the status of a job.

    :param jobid: the jobid returned by 'submit'
    :return: Either PENDING, STARTED, SUCCESS, FAILURE, RETRY, or REVOKED.
    """

    _log.info ("endpoints.status request for job %s" % jobid)

    from hommod_rest.application import celery
    result = celery.AsyncResult(jobid)
    job_status = result.status

    response = {'status': job_status}
    if result.failed():
        response ['message'] = str (result.traceback)

    _log.debug ("endpoints.status: response for job %s: %s" %
                (jobid, str (job_status)))

    return jsonify(response)


@bp.route('/result/<job_id>/', methods=['GET'])
def result(job_id):
    from hommod_rest.application import celery
    result = celery.AsyncResult(job_id)
    path = result.result
    if not path:
        return jsonify({'model_created': False})
    else:
        return jsonify({'model_created': True})


@bp.route('/get_model_file/<jobid>.pdb', methods=['GET'])
@bp.route('/get_model_file/<jobid>.PDB', methods=['GET'])
def get_model_file(jobid):

    """
    Get the pdb file, created by the modeling job.

    :param jobid: the jobid returned by 'submit'
    :return: The pdb file created by the job. If the job status is not SUCCESS, this method returns an error.
    """

    _log.info("endpoints: model request for job %s" % jobid)

    from hommod_rest.application import celery
    result = celery.AsyncResult(jobid)
    path = result.result
    if not path:
        # no model could be created
        message = 'no model was created for job %s' % jobid
        _log.warn (message)
        return jsonify({'error': message}), 400

    try:
        contents = extract_model(path)
    except Exception as e:
        error = 'failed to get all data from %s: %s' % (path, str (e))
        _log.error(error)
        return jsonify({'error': error}), 500

    _log.debug("endpoints: model successfully retrieved for job %s" % jobid)

    return Response(contents, mimetype='chemical/x-pdb')


@bp.route('/get_metadata/<jobid>/', methods=['GET'])
def get_metadata(jobid):

    """
    Get the metadata of the model, created by the modeling job.

    :param jobid: the jobid returned by 'submit'
    :return: The json metadata object. If the job status is not SUCCESS,
             this method returns an error.
    """

    _log.debug("metadata request for job %s" % jobid)

    from hommod_rest.application import celery
    result = celery.AsyncResult(jobid)
    path = result.result
    if not path:
        return jsonify({})

    try:
        data = extract_info(path)
        data['alignment'] = extract_alignment(path)
    except Exception as e:
        error = 'failed to get all data from %s: %s' % (path, str (e))
        _log.error(error)
        return jsonify({'error': error}), 500

    _log.debug("endpoints: metadata successfully retrieved for job %s" % jobid)

    return jsonify(data)
