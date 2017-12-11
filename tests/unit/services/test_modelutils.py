from nose.tools import eq_

from hommod_rest.services.modelutils import getChainCAsSeqSecStr


def test_getChainCAsSeqSecStr():

    class FakeYasaraModule(object):
        def ListAtom(self, s, f):
            return ["3524 XPC", "3527 ALA", "3537 XCP", "3541 ARG"]
        def SecStrRes(self, s):
            return ["H", "H", "H", "H"]

    module = FakeYasaraModule()

    CAs, seq, secstr = getChainCAsSeqSecStr(module, 1, 'D')
    eq_(seq, "XAXR")
