import json
import logging
from time import sleep

from mock import patch
from nose.tools import with_setup, ok_, eq_

from hommod.application import app as flask_app


_log = logging.getLogger(__name__)


def setup():
    flask_app.testing = True
    global flask_client
    flask_client = flask_app.test_client()

def teardown():
    pass


@patch('hommod.tasks.create_model.apply_async')
@patch('hommod.application.celery.AsyncResult')
@patch('hommod.controllers.storage.model_storage.extract_model')
@patch('hommod.controllers.storage.model_storage.extract_alignments')
@patch('hommod.controllers.storage.model_storage.extract_selected_targets')
@with_setup(setup, teardown)
def test_interface(mock_targets, mock_alignments, mock_model, mock_async, mock_result):

    target_id = 'target'

    class FakeResult:
        def __init__(self, job_id):
            self.status = 'SUCCESS'
            self.task_id = str(job_id)

        def get(self):
            return [target_id]

        def failed(self):
            return False

    mock_result.return_value = FakeResult('no-job')
    mock_async.return_value = FakeResult('no-job')

    mock_alignments.return_value = []
    mock_targets.return_value = {'A': target_id}
    mock_model.return_value = ''

    sequence = "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"

    r = flask_client.post('/api/submit/',
                          data={'sequence': sequence,
                                'position': 4,
                                'species_id': 'CRAAB',
                                'template_id': '1CRN_A'})
    eq_(r.status_code, 200)

    job_id = json.loads(r.data)['jobid']

    while True:
        r = flask_client.get('/api/status/%s/' % job_id)
        eq_(r.status_code, 200)

        status = json.loads(r.data)['status']

        _log.debug("got status {}".format(status))

        ok_(status not in ['FAILURE', 'REVOKED'])

        if status == 'SUCCESS':
            break
        else:
            sleep(1)

    _log.debug("getting y/n")

    r = flask_client.get('/api/result/%s/' % job_id)
    eq_(r.status_code, 200)
    ok_(json.loads(r.data)['model_created'])

    _log.debug("getting model")

    r = flask_client.get('/api/get_model_file/%s.pdb' % job_id)
    eq_(r.status_code, 200)

    _log.debug("getting metadata")

    r = flask_client.get('/api/get_metadata/%s/' % job_id)
    eq_(r.status_code, 200)
