from hommod_rest.services.modelutils import get_pdb_contents
from nose.tools import ok_

def test_get_pdb_contents():
    s = get_pdb_contents('1crn')
    ok_(len(s) > 0)
    # must not be binary:
    ok_(all(ord(c) < 128 for c in s))
