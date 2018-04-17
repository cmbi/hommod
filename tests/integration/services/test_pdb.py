from nose.tools import ok_

from hommod.services.pdb import get_pdb_contents


def test_get_pdb_contents():
    contents = get_pdb_contents('5gox')

    ok_(len(contents) > 0)
