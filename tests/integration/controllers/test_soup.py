import os
import tempfile
import shutil
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
            work_dir = tempfile.mkdtemp()

            try:
                soup.yasara.CD(work_dir)
                soup.yasara.LoadPDB('1crn', download='yes')
                eq_(soup.yasara.ListMol('all', 'MOL'), ['A'])

                soup.yasara.Exit()
            finally:
                shutil.rmtree(work_dir)

    ps = [TestProcess(), TestProcess(), TestProcess()]
    for p in ps:
        p.start()
    for p in ps:
        p.join()
        eq_(p.exitcode, 0)


@with_setup(setup, teardown)
def test_yasara_restart():
    work_dir = tempfile.mkdtemp()

    try:
        soup.yasara.CD(work_dir)
        soup.yasara.LoadPDB('1crn', download='yes')
        soup.yasara.Exit()

        soup.yasara.CD(work_dir)
        soup.yasara.LoadPDB('1crn', download='yes')
        eq_(soup.yasara.ListMol('all', 'MOL'), ['A'])
    finally:
        shutil.rmtree(work_dir)
