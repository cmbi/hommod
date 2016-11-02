from nose.tools import with_setup
from hommod_rest.services.model import modeler
import hommod_rest.default_settings as config


def _setup():
    modeler.yasara_dir = config.YASARADIR
    modeler.execution_root_dir = config.EXECUTIONDIR
    modeler.model_root_dir = config.MODELDIR
    modeler.template_blacklist = config.TEMPLATE_BLACKLIST


def _tear_down():
    pass


@with_setup(_setup, _tear_down)
def test_set_template():
    obj, oligomerization = modeler._set_template("3DS1")

    assert obj is not None
    assert oligomerization >= 1
