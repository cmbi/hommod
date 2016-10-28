from hommod_rest.factory import create_app, create_celery_app

import logging
_log = logging.getLogger (__name__)

def test_create_model():

    app = create_app ({'TESTING': True, 'CELERY_ALWAYS_EAGER': True})
    celery = create_celery_app (app)

    # Best pick a simple model which doesn't take much time ..
    from hommod_rest.tasks import create_model
    result = create_model.delay (
                "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN",
                "CRAAB", 25)

    path = result.get ()
    _log.debug ("Recieved model path: \'%s\'" % path)

    assert (len (path) > 0)
