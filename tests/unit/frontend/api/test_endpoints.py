import json

from celery.result import AsyncResult
from mock import create_autospec, Mock, patch
from nose.tools import eq_

from hommod_rest.factory import create_app


class TestEndpoints:
    def setup(self):
        self.app = create_app({'TESTING': True}).test_client()

    def teardown(self):
        pass

    def test_submit_all_none(self):
        data = {
            'sequence': None, 'position': None,
            'species_id': None, 'template_id': None,
        }
        rv = self.app.post('/api/submit/', data=data)

        eq_(rv.status_code, 400)
        eq_(json.loads(rv.data)['error'], 'invalid input')

    @patch('hommod_rest.tasks.create_model')
    def test_submit_template_none(self, mock_create_model):
        mock_result = create_autospec(AsyncResult)
        mock_result.status = 'PENDING'
        mock_result.task_id = '1234567890'
        mock_create_model.apply_async.return_value = mock_result

        data = {
            'sequence': 'AMELK', 'position': 10,
            'species_id': 'HUMAN', 'template_id': None,
        }
        rv = self.app.post('/api/submit/', data=data)

        eq_(json.loads(rv.data)['jobid'], '1234567890')
