import logging

from flask import Blueprint, jsonify, request, Response

from hommod_rest.managejobs import create_model
from hommod_rest.services.utils import (extract_info, extract_alignment,
                                        extract_model)

_log = logging.getLogger(__name__)


bp = Blueprint('hommod', __name__, url_prefix='/api')


@bp.route('/submit/', methods=['POST'])
def submit():
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

    jobid = create_model(sequence, species_id, position)
    print jobid

    return jsonify({'jobid': jobid})


@bp.route('/status/<jobid>/', methods=['GET'])
def status(jobid):
    from hommod_rest.application import celery
    result = celery.AsyncResult(jobid)
    return jsonify({'status': result.status})


@bp.route('/get_model_file/<jobid>.pdb', methods=['GET'])
def get_model_file(jobid):

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
def get_metadata(jobid):

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
