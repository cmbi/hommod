from hommod_rest.services.align import aligner

from nose.tools import ok_, with_setup
import hommod_rest.default_settings as config

def _setup ():

    aligner.clustalExe = config.CLUSTAL
    aligner.kmadExe = config.KMAD

def _tear_down ():
    pass

@with_setup(_setup, _tear_down)
def test_clustal ():

    toalign = {"s1": "AAACCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGSDYAN",
               "s2": "TTCCPSIVASNVCRLPGTPEAICATYTGCIIIPGATCPGDYA"}

    aligned = aligner.clustalAlign (toalign)
    assert "s1" in aligned
    assert "s2" in aligned

@with_setup(_setup, _tear_down)
def test_kmad ():

    targetseq = "AAACCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGSDYAN"
    templateseq    = "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"
    templatesecstr = " EE SSHHHHHHHHHHHTTT  HHHHHHHHS EE SSS   GGG  "

    aligned = aligner.kmadAlign (templateseq, templatesecstr, targetseq)
    assert "target" in aligned
    assert "template" in aligned
