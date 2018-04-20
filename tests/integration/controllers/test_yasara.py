import os
import sys
import tempfile
import shutil
from multiprocessing import Process
from threading import Thread
import logging

from mock import patch
from nose.tools import ok_, eq_, with_setup

from hommod.default_settings import YASARA_DIR
from hommod.controllers.yasara import YasaraContext


_log = logging.getLogger(__name__)

def setup():
    pass


def teardown():
    pass


@with_setup(setup, teardown)
def test_yasara_multithread():
    class TestThread(Thread):
        def __init__(self):
            Thread.__init__(self)
            self._ex = None

        def run(self):
            work_dir = tempfile.mkdtemp()
            yasara = YasaraContext(YASARA_DIR)

            try:
                yasara.CD(work_dir)
                yasara.DownloadPDB('1crn')
                eq_(yasara.ListMol('all', 'MOL'), ['A'])
            except:
                self._ex = sys.exc_info()
            finally:
                shutil.rmtree(work_dir)

        def join(self):
            Thread.join(self)

            if self._ex is not None:
                raise self._ex[0], self._ex[1], self._ex[2]


    ts = [TestThread(), TestThread(), TestThread()]
    for t in ts:
        t.start()
    for t in ts:
        t.join()


@with_setup(setup, teardown)
def test_yasara_multiprocess():
    class TestProcess(Process):
        def run(self):
            work_dir = tempfile.mkdtemp()
            yasara = YasaraContext(YASARA_DIR)

            try:
                yasara.CD(work_dir)
                yasara.DownloadPDB('1crn')
                eq_(yasara.ListMol('all', 'MOL'), ['A'])
            finally:
                shutil.rmtree(work_dir)

    ps = [TestProcess(), TestProcess(), TestProcess()]
    for p in ps:
        p.start()
    for p in ps:
        p.join()
        eq_(p.exitcode, 0)


@patch('hommod.controllers.yasara.YasaraContext.__del__')
@with_setup(setup, teardown)
def test_yasara_shutdown(mock_del):

    work_dir = tempfile.mkdtemp()
    os.chdir(work_dir)
    yasara = YasaraContext(YASARA_DIR)

    try:
        yasara.CD(work_dir)
        del yasara

        ok_(mock_del.called)
        ok_(not os.path.isfile(os.path.join(work_dir, 'errorexit.txt')))
    finally:
        shutil.rmtree(work_dir)
