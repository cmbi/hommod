from nose.tools import with_setup
from mock import patch

from hommod_rest.services.interpro import interpro
from hommod_rest.default_settings import INTERPRO_TIMEOUT


def setup():
    interpro.storage_dir = '/does/not/exist'
    interpro.job_timeout = INTERPRO_TIMEOUT

def teardown():
    pass


@with_setup(setup, teardown)
@patch("os.path.isdir")
@patch("filelock.FileLock.acquire")
@patch("bz2.BZ2File")
@patch("hommod_rest.services.interpro._interpro_submit")
@patch("hommod_rest.services.interpro._interpro_get_status")
@patch("hommod_rest.services.interpro._interpro_get_result")
def test_interpro(mock_result, mock_status, mock_submit,
                  mock_bz2, mock_lock_acquire, mock_isdir):

    def fake_lock():
        pass

    jobid = 'hi'
    status = 'FINISHED'
    result = '<protein-matches><protein></protein></protein-matches>'

    class FakeBZ2(object):
        def read(self):
            return result
        def write(self, data):
            pass

    mock_bz2.return_value = FakeBZ2()
    mock_result.return_value = result
    mock_status.return_value = status
    mock_submit.return_value = jobid
    mock_isdir.return_value = True
    mock_lock_acquire.side_effect = fake_lock

    interpro.get_domain_locations("HI")
