from hommod_rest.services.domainalign import getAlignments
from hommod_rest.services.modelutils import getInterproDomainLocations
from hommod_rest.factory import create_app, create_celery_app
from nose.tools import ok_, with_setup

from hommod_rest.services.interpro import interpro


def setupApps():
    app = create_app({'TESTING': True, 'DEBUG':True, \
                      'CELERY_ALWAYS_EAGER': True, 'RETRY_FAILURE' : True,
                      'YASARADIR': '/home/cbaakman/Yasara/yasara'})
    celery_app = create_celery_app(app)


def tearDown():
    pass


@with_setup(setupApps,tearDown)
def test_COX41_BOVIN():
    seq="""
    MLATRVFSLIGRRAISTSVCVRAHGSVVKSEDYALPSYVDRRDYPLPDVAHVKNLSASQK
    ALKEKEKASWSSLSIDEKVELYRLKFKESFAEMNRSTNEWKTVVGAAMFFIGFTALLLIW
    EKHYVYGPIPHTFEEEWVAKQTKRMLDMKVAPIQGFSAKWDYDKNEWKK
    """.replace("\n",'')

    domainranges=interpro.getInterproDomainLocations(seq)

    alignments = getAlignments(domainranges, seq)

    ok_(len(alignments) > 0)
