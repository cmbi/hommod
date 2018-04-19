from multiprocessing import Process
import logging

from nose.tools import ok_, eq_, with_setup

from hommod.default_settings import YASARA_DIR
from hommod.controllers.soup import soup


_log = logging.getLogger(__name__)

def setup():
    soup.yasara_dir = YASARA_DIR

def teardown():
    pass


@with_setup(setup, teardown)
def test_yasara_multiprocess():
    class TestProcess(Process):
        def run(self):
            soup.yasara.LoadPDB('1crn', download='yes')
            ok_(len(soup.yasara.ListMol('all')) > 0)

            soup.yasara.Exit()

    ps = [TestProcess(), TestProcess(), TestProcess()]
    for p in ps:
        p.start()
    for p in ps:
        p.join()
        eq_(p.exitcode, 0)
