from hommod_rest.services.align import aligner

from nose.tools import with_setup
import hommod_rest.default_settings as config


def setup():
    aligner.clustal_exe = config.CLUSTAL
    aligner.kmad_exe = config.KMAD


def teardown():
    pass


@with_setup(setup, teardown)
def test_clustal():

    toalign = {"s1": "AAACCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGSDYAN",
               "s2": "TTCCPSIVASNVCRLPGTPEAICATYTGCIIIPGATCPGDYA"}

    aligned = aligner.clustal_align(toalign)
    assert "s1" in aligned
    assert "s2" in aligned


@with_setup(setup, teardown)
def test_kmad():

    targetseq = "AAACCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGSDYAN"
    templateseq = "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"
    templatesecstr = " EE SSHHHHHHHHHHHTTT  HHHHHHHHS EE SSS   GGG  "

    aligned = aligner.kmad_align(templateseq, templatesecstr, targetseq)
    assert "target" in aligned
    assert "template" in aligned
