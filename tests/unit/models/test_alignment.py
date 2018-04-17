from nose.tools import eq_

from hommod.models.align import BlastAlignment, TargetTemplateAlignment, Alignment


def test_get_hit_accession_code():
    alignment = BlastAlignment('sp|P00395|COX1_HUMAN', '', '', 0, 0, '', 0, 0, '')
    eq_(alignment.get_hit_accession_code(), 'P00395')


def test_get_hit_chain_id():
    alignment = BlastAlignment('pdb|1crn|A', '', '', 0, 0, '', 0, 0, '')
    eq_(alignment.get_hit_chain_id(), 'A')


def test_get_hit_type():
    alignment = BlastAlignment('pdb|1crn|A', '', '', 0, 0, '', 0, 0, '')
    eq_(alignment.get_hit_type(), 'pdb')

    alignment = BlastAlignment('sp|P00395|COX1_HUMAN', '', '', 0, 0, '', 0, 0, '')
    eq_(alignment.get_hit_type(), 'sp')

    alignment = BlastAlignment('1crn_A', '', '', 0, 0, '', 0, 0, '')
    eq_(alignment.get_hit_type(), 'pdb')


def test_get_percentage_coverage():
    alignment = BlastAlignment('pdb|1crn|A', 'AVAVAVAVAV', '',
                               1, 10, 'AVAVAVAVAV',
                               1, 10, 'A-A-A-A-A-')
    eq_(alignment.get_percentage_coverage(), 50.0)

    alignment = TargetTemplateAlignment('AVAVAVAVAV',
                                        'A-A-A-A-A-')
    eq_(alignment.get_percentage_coverage(), 50.0)


def test_get_percentage_identity():
    alignment = BlastAlignment('pdb|1crn|A', 'AVAVAVAVAV', '',
                               1, 10, 'AVAVAVAVAV',
                               1, 10, 'ATATATATAT')
    eq_(alignment.get_percentage_identity(), 50.0)

    alignment = TargetTemplateAlignment('AVAVAVAVAV',
                                        'ATATATATAT')
    eq_(alignment.get_percentage_identity(), 50.0)


def test_count_aligned_residues():
    alignment = BlastAlignment('pdb|1crn|A', 'AVAVAVAVAV', '',
                               1, 10, 'AVAVAVAVAV',
                               1, 10, 'A-A-A-A-A-')
    eq_(alignment.count_aligned_residues(), 5)

    alignment = TargetTemplateAlignment('AVAVAVAVAV',
                                        'A-A-A-A-A-')
    eq_(alignment.count_aligned_residues(), 5)


def test_get_covered_template_residues_indices():
    alignment = TargetTemplateAlignment('---VA-AVAV',
                                        'A-AAAAA-A-')

    indices = alignment.get_covered_template_residues_indices()
    eq_(indices, [2, 3, 5, 6])


def test_get_relative_span():
    alignment = TargetTemplateAlignment('---VA-AVAV-',
                                        '--AAAAA-A-A')
    r = alignment.get_relative_span()
    eq_(r.start, 1)
    eq_(r.end, 6)


def test_get_query_range():
    alignment = BlastAlignment('pdb|1crn|A', 'AVAVAVAVAV', '',
                               3, 7, 'AVAVA',
                               4, 8, 'ATATA')
    r = alignment.get_query_range()
    eq_(r.start, 2)
    eq_(r.end, 7)


def test_as_dict():
    d = {'1': "AVAVAVAVA",
         '2': "ATATATATA"}
    alignment = Alignment(d)
    r = alignment.as_dict()

    eq_(r, d)
