from hommod_rest.services.domainalign import getAlignments
from hommod_rest.services.modelutils import getInterproDomainLocations
from hommod_rest.factory import create_app, create_celery_app
from nose.tools import ok_, with_setup

from hommod_rest.services.modelproc import modeler, findOrthologsToSeq

import logging
_log = logging.getLogger(__name__)


def setupApps():
    app = create_app({'TESTING': True, 'DEBUG':True, \
                      'CELERY_ALWAYS_EAGER': True, 'RETRY_FAILURE' : True,
                      'YASARADIR': '/home/cbaakman/Yasara/yasara'})
    celery_app = create_celery_app(app)


def tearDown():
    pass


@with_setup(setupApps,tearDown)
def test_1OCO_BOVIN():
    tempobj,oligomerisation=modeler.setTemplate('1OCO')
    chainOrder,templateChainSequences=modeler.getChainOrderAndSeqs(tempobj)

    for chain in chainOrder:
        _log.debug("finding orthologs for %s:\n%s" % (chain,templateChainSequences[chain]) )
        orthologs = findOrthologsToSeq(templateChainSequences[chain],'BOVIN')

        ok_( len(orthologs) > 0 )

