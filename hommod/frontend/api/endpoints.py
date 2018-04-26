import inspect
import logging
import re
import os
import traceback

from flask import Blueprint, render_template, request, jsonify, Response

from hommod.models.template import TemplateID
from hommod.controllers.storage import model_storage
from hommod.controllers.sequence import is_protein_sequence
from hommod.controllers.method import select_best_model

bp = Blueprint('api', __name__, url_prefix='/api')

_log = logging.getLogger(__name__)


_P_TEMPLATE_ID = re.compile(r'[0-9][0-9a-zA-Z]{3}_[0-9a-zA-Z]')

def _validate_input_data(form):
    if 'sequence' not in form:
        raise ValueError("Missing sequence input")
    elif not is_protein_sequence(form['sequence']):
        raise ValueError("Invalid sequence data")
    sequence = form['sequence'].upper()

    if 'species_id' not in form:
        raise ValueError("Missing species_id input")
    elif not form['species_id'].isalpha():
        raise ValueError("Invalid species_id data")
    species_id = form['species_id']

    if 'position' in form:
        position = int(form['position'])

        if position < 1 or position > len(sequence):
            raise ValueError("Residue position out of range")
    else:
        position = None

    if 'template_id' in form:
        if not _P_TEMPLATE_ID.match(form['template_id']):
            raise ValueError("Invalid template id data")
        pdbid, chain_id = form['template_id'].split('_')
        template_id = TemplateID(pdbid, chain_id)
    else:
        template_id = None

    return sequence, species_id, position, template_id


@bp.route('/submit/', methods=['POST'])
def submit():

    """ 
    Request a model for the given parameters.

    :param sequence: target sequence for the model
    :param species_id: uniprot species id for the model
    :param position: optional position of the required residue in the sequence, starting 1
    :param template_id: optional pdbid and chain id, separated by '_'
    :return: a json object, containing the field 'jobid' or an error if the input is incorrect
    """

    try:
        sequence, species_id, position, template_id = _validate_input_data(request.form)
    except Exception as e:
        return jsonify({'error': str(e)}), 400

    from hommod.tasks import create_model
    result = create_model.apply_async((sequence, species_id, position, template_id))

    return jsonify({'jobid': result.task_id})


@bp.route('/get_model_if_exists/', methods=['POST'])
def get_model_if_exists():

    """
    Get a model if it exists in the set of earlier created models.

    :param sequence: target sequence for the model
    :param species_id: uniprot species id for the model
    :param position: optional position of the required residue in the sequence, starting 1
    :param template_id: optional 'pdbid'-'chain' for the template
    :return: a json object, containing the field 'model_ids'
    """

    try:
        sequence, species_id, position, template_id = _validate_input_data(request.form)
    except Exception as e:
        return jsonify({'error': str(e)}), 400

    paths = model_storage.list_models(sequence, species_id, position, template_id)

    model_ids = [model_storage.get_model_name_from_path(path) for path in paths]

    return jsonify({'model_ids': model_ids})


@bp.route('/status/<job_id>/', methods=['GET'])
def status(job_id):

    """
    Request the status of a job.

    :param jobid: the job id returned by 'submit'
    :return: Either PENDING, STARTED, SUCCESS, FAILURE, RETRY, or REVOKED.
    """

    from hommod.application import celery
    result = celery.AsyncResult(job_id)

    response = {'status': result.status}
    if result.failed():
        response['message'] = str(result.traceback)

    return jsonify(response)


@bp.route('/result/<job_id>/', methods=['GET'])
def result(job_id):

    """
    Request whether a job has created a model or not.

    :param jobid: the job id returned by 'submit'
    :return: a json object containing the boolean field 'model_created'
    """

    from hommod.application import celery
    result = celery.AsyncResult(job_id)

    if result.status != 'SUCCESS':
        return jsonify({'error': "{} has status {}".format(job_id, result.status)}), 400

    try:
        path = result.get()
    except:
        return jsonify({'model_created': False})

    if path is not None:
        return jsonify({'model_created': True})
    else:
        return jsonify({'model_created': False})


@bp.route('/get_model_file/<job_id>.pdb', methods=['GET'])
@bp.route('/get_model_file/<job_id>.PDB', methods=['GET'])
def get_model_file(job_id):

    """ 
    Get the pdb file, created by the modeling job.

    :param jobid: the job_id returned by 'submit'
    :return: The pdb file created by the job. If the job status is not SUCCESS,
             this method returns an error.
    """

    from hommod.application import celery
    result = celery.AsyncResult(job_id)

    if result.status != 'SUCCESS':
        return jsonify({'error': "{} has status {}".format(job_id, result.status)}), 400

    try:
        path = result.get()
    except Exception:
        return jsonify({'error': result.traceback}), 500

    if path is None:
        message = 'Job %s finished, but without creating a model. This could be due to lack of a suitable template.' % job_id
        return jsonify({'error': message}), 500

    try:
        contents = model_storage.extract_model(path)
        return Response(contents, mimetype='chemical/x-pdb')
    except:
        return jsonify({'error': traceback.format_exc()}), 500


@bp.route('/get_model_file_by_model_id/<model_id>.pdb', methods=['GET'])
@bp.route('/get_model_file_by_model_id/<model_id>.PDB', methods=['GET'])
def get_model_file_by_model_id(model_id):

    """
    Get the pdb file, created by the modeling job.

    :param model_id: the id returned by 'get_model_if_exists'
    :return: The pdb file created by the job.
    """

    path = model_storage.get_tar_path_from_name(model_id)
    if not os.path.isfile(path):
        return jsonify({'error': "no such model"}), 400

    try:
        contents = model_storage.extract_model(path)
        return Response(contents, mimetype='chemical/x-pdb')
    except:
        message = str(e)
        return jsonify({'error': traceback.format_exc()}), 500


@bp.route('/get_metadata/<job_id>/', methods=['GET'])
def get_metadata(job_id):

    """ 
    Get the metadata of the model, created by the modeling job.

    :param jobid: the job_id returned by 'submit'
    :return: The json metadata object. If the job status is not SUCCESS,
             this method returns an error.
    """

    from hommod.application import celery
    result = celery.AsyncResult(job_id)

    if result.status != 'SUCCESS':
        return jsonify({'error': "{} has status {}".format(job_id, result.status)}), 400

    try:
        path = result.get()
    except Exception:
        return jsonify({'error': result.traceback}), 500

    if path is None:
        message = 'Job %s finished, but without creating a model. This could be due to lack of a suitable template.' % job_id
        return jsonify({'error': message}), 500

    try:
        data = {}
        data['selected_targets'] = model_storage.extract_selected_targets(path)
        data['alignments'] = [alignment.as_dict()
                              for alignment in model_storage.extract_alignments(path)]

        return jsonify(data)
    except:
        return jsonify({'error': traceback.format_exc()}), 500


@bp.route('/get_metadata_by_model_id/<model_id>/', methods=['GET'])
def get_metadata_by_model_id(model_id):

    """
    Get the metadata of the model, created by the modeling job.

    :param model_id: the id returned by 'get_model_if_exists'
    :return: The json metadata object.
    """

    path = model_storage.get_tar_path_from_name(model_id)
    if not os.path.isfile(path):
        return jsonify({'error': "no such model"}), 400

    try:
        data = {}
        data['selected_targets'] = model_storage.extract_selected_targets(path)
        data['alignments'] = [alignment.as_dict()
                              for alignment in model_storage.extract_alignments(path)]

        return jsonify(data)
    except:
        return jsonify({'error': traceback.format_exc()}), 500


@bp.route('/', methods=['GET'])
def api_doc():
    fs = [submit, status, result, get_model_file, get_metadata,
          get_model_if_exists, get_model_file_by_model_id, get_metadata_by_model_id]
    docs = []
    for f in fs:
        src = inspect.getsourcelines(f)
        m = re.search(r"@bp\.route\(\'(.*?)\'\, methods\=\['([A-Z]*)']\)",
                      src[0][0])
        if not m:  # pragma: no cover
            _log.debug("Unable to document function '{}'".format(f))
            continue

        url = m.group(1)
        method = m.group(2)
        docstring = inspect.getdoc(f)
        docs.append((url, method, docstring))
    return render_template('api/docs.html', docs=docs)


@bp.errorhandler(Exception)
def exception_error_handler(error):  # pragma: no cover
    tb = traceback.format_exc()
    _log.error("Unhandled exception:\n{}".format(tb))
    return jsonify({'error': tb}), 500
