from nose.tools import eq_

from hommod.models.align import BlastAlignment


def test_get_hit_accession_code():
    alignment = BlastAlignment('sp|P00395|COX1_HUMAN', 0, 0, '', 0, 0, '')
    eq_(alignment.get_hit_accession_code(), 'P00395')
