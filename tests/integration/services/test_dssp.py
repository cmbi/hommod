from nose.tools import with_setup, ok_

from hommod.services.dssp import dssp
from hommod.models.template import TemplateID
from hommod import default_settings as settings


def setup():
    dssp.dssp_dir = settings.DSSP_DIR

def teardown():
    pass


@with_setup(setup, teardown)
def test_3jb9():
    ok_(dssp.has_secondary_structure(TemplateID('3jb9', 'B')))
