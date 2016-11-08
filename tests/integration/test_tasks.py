from hommod_rest.factory import create_app, create_celery_app

import logging
_log = logging.getLogger(__name__)


class TestTasks:

    @classmethod
    def setup_class(self):
        self.app = create_app({'TESTING': True, 'CELERY_ALWAYS_EAGER': True})
        self.celery = create_celery_app(self.app)

    def test_create_model(self):
        # Best pick a simple model which doesn't take much time ..
        from hommod_rest.tasks import create_model
        path = create_model("IICCASITARSDFDVCRLAGTAQAVCAVFVGCVIVPAATCPGDFGD",
                            "CRAAB", 25, None)

        _log.debug("Recieved model path: \'%s\'" % path)

        assert (len(path) > 0)

    def test_take_template(self):
        # Best pick a simple model which doesn't take much time ..
        from hommod_rest.tasks import create_model
        path = create_model("TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN",
                            "CRAAB", 25, None)

        _log.debug("Recieved model path: \'%s\'" % path)

        assert (len(path) > 0)
