from nose.tools import with_setup, ok_

from hommod.controllers.clustal import clustal_aligner
from hommod.default_settings import CLUSTALW_EXE


def setup():
    clustal_aligner.clustalw_exe = CLUSTALW_EXE


def teardown():
    pass


@with_setup(setup, teardown)
def test_clustal():
    toalign = {"s1": "AAACCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGSDYAN",
               "s2": "TTCCPSIVASNVCRLPGTPEAICATYTGCIIIPGATCPGDYA"}

    aligned = clustal_aligner.align(toalign)
    ok_("s1" in aligned.aligned_sequences)
    ok_("s2" in aligned.aligned_sequences)

