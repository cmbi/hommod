from nose.tools import with_setup
from hommod_rest.services.model import modeler
import hommod_rest.default_settings as config


def setup():
    modeler.yasara_dir = config.YASARADIR
    modeler.execution_root_dir = config.EXECUTIONDIR
    modeler.model_root_dir = config.MODELDIR
    modeler.template_blacklist = config.TEMPLATE_BLACKLIST


def teardown():
    pass


@with_setup(setup, teardown)
def test_set_template():
    obj, oligomerization = modeler._set_template("3DS1")

    assert obj is not None
    assert oligomerization >= 1
