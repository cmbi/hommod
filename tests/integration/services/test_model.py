from nose.tools import eq_

from hommod_rest.services.model import modeler
from hommod_rest.default_settings import YASARADIR, TEMPLATE_BLACKLIST
modeler.yasara_dir = YASARADIR
modeler.execution_root_dir = '/tmp'
modeler.model_root_dir = '/tmp'
modeler.template_blacklist = TEMPLATE_BLACKLIST


def test_fix_chain_breaks_5MQF():
    modeler._set_template('5MQF')
    eq_(1, len(modeler.yasara.ListMol('l')))


