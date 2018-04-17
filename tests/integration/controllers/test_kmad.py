from nose.tools import with_setup, eq_, ok_

from hommod.controllers.kmad import kmad_aligner
from hommod.default_settings import KMAD_EXE


def setup():
    kmad_aligner.kmad_exe = KMAD_EXE


def teardown():
    pass


@with_setup(setup, teardown)
def test_kmad():
    target_seq = "AAACCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGSDYAN"
    template_seq = "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"
    template_secstr = " EE SSHHHHHHHHHHHTTT  HHHHHHHHS EE SSS   GGG  "

    aligned = kmad_aligner.align(template_seq, template_secstr, target_seq)
    ok_(len(aligned.target_alignment) > 0)
    eq_(len(aligned.target_alignment), len(aligned.template_alignment))


@with_setup(setup, teardown)
def test_kmad_X():
    target_seq = "AAAAAAAAAAAAAAA"
    template_seq = "XAXRXLXKXGDAFNR"
    template_secstr = "               "

    aligned = kmad_aligner.align(template_seq, template_secstr, target_seq)
    eq_(aligned.template_alignment.replace('-',''), template_seq)
