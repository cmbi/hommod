from nose.tools import ok_, with_setup

from hommod_rest.services.model import modeler, findOrthologsToSeq
from hommod_rest.services.blast import blaster
import hommod_rest.default_settings as config

import logging
_log = logging.getLogger(__name__)


def setupApps():

    modeler.yasara_dir = config.YASARADIR
    modeler.execution_root_dir = config.EXECUTIONDIR
    modeler.model_root_dir = config.MODELDIR
    modeler.template_blacklist = config.TEMPLATE_BLACKLIST

    blaster.uniprotDB = config.UNIPROTDB
    blaster.blastpExe = config.BLASTP


def tearDown():
    pass


@with_setup(setupApps, tearDown)
def test_1OCO_BOVIN():
    tempobj, oligomerisation = modeler._set_template('1OCO')
    chainOrder, templateChainSequences = modeler.getChainOrderAndSeqs(tempobj)

    for chain in chainOrder:
        _log.debug("finding orthologs for %s:\n%s"
                   % (chain, templateChainSequences[chain]))
        orthologs = findOrthologsToSeq(templateChainSequences[chain], 'BOVIN')

        ok_(len(orthologs) > 0)
