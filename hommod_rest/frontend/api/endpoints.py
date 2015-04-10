import logging

from flask import Blueprint, jsonify, request

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


@bp.route('/results/<jobid>/', methods=['GET'])
def results(jobid):
    from hommod_rest.application import celery
    result = celery.AsyncResult(jobid)

    if result.ready():

        paths = result.result
        results = []
        for path in paths:

            try:
                info = extract_info(path)
                alignment = extract_alignment(path)
                contents = extract_model(path)
            except:
                _log.warn('failed to get all data from %s' % path)
                continue

            results.append({
                'model': contents, 'alignment': alignment, 'info': info})

        return jsonify({'results': results})
    else:
        return jsonify({'error': 'not ready yet'}), 400
