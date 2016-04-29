from hommod_rest.services.domainalign import getAlignments
from nose.tools import ok_, with_setup

from hommod_rest.services.interpro import interpro
from hommod_rest.services.blast import blaster
from hommod_rest.services.secstr import secstr

import hommod_rest.default_settings as config


def setupApps():
    interpro.interproExe = config.INTERPROSCAN
    interpro.storageDir = config.INTERPRODIR
    blaster.uniprotDB = config.UNIPROTDB
    blaster.templatesDB = config.TEMPLATESDB
    blaster.blastpExe = config.BLASTP
    secstr.dssp_dir = config.DSSPDIR
    secstr.yasara_dir = config.YASARADIR


def tearDown():
    pass


@with_setup(setupApps, tearDown)
def test_COX41_BOVIN():
    seq = """
    MLATRVFSLIGRRAISTSVCVRAHGSVVKSEDYALPSYVDRRDYPLPDVAHVKNLSASQK
    ALKEKEKASWSSLSIDEKVELYRLKFKESFAEMNRSTNEWKTVVGAAMFFIGFTALLLIW
    EKHYVYGPIPHTFEEEWVAKQTKRMLDMKVAPIQGFSAKWDYDKNEWKK
    """.replace("\n", '')

    domainranges = interpro.getInterproDomainLocations(seq)

    alignments = getAlignments(domainranges, seq)

    ok_(len(alignments) > 0)
