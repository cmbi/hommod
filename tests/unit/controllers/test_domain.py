from nose.tools import eq_, ok_, with_setup

from hommod.controllers.domain import domain_aligner
from hommod.models.range import SequenceRange
from hommod.models.align import TargetTemplateAlignment, BlastAlignment
from hommod.models.template import TemplateID


def setup():
    from hommod import default_settings as settings
    domain_aligner.similar_ranges_min_overlap_percentage = settings.SIMILAR_RANGES_MIN_OVERLAP_PERCENTAGE
    domain_aligner.similar_ranges_max_length_difference_percentage = settings.SIMILAR_RANGES_MAX_LENGTH_DIFFERENCE_PERCENTAGE
    domain_aligner.forbidden_interpro_domains = settings.FORBIDDEN_INTERPRO_DOMAINS
    domain_aligner.highly_homologous_percentage_identity = settings.HIGHLY_HOMOLOGOUS_PERCENTAGE_IDENTITY
    domain_aligner.template_blast_databank = settings.TEMPLATE_BLAST_DATABANK
    domain_aligner.min_percentage_coverage = settings.DOMAIN_MIN_PERCENTAGE_COVERAGE


def teardown():
    pass

def test_get_template_sequence_in_target_range():

    alignment = TargetTemplateAlignment(
       "PHTSHSWLCDGRLLCLHDPSNKNNWKIFRECWKQGQPVLVSGVHKKLK" +
       "SELWKPEAFSQEFGDQDVDLVNCRNCAIISDVKVRDFWDGFEIICKRL" +
       "RSEDGQPMVLKLKDWPPGEDFRDMMPTRFEDLMENLPLPEYTKRDGRL" +
       "NLASRLPSYFVRPDLGPKMYNAYGLITAEDRRVGTTNLHLDVSDAVNV" +
       "MVYVGIPIGEG-AHDEEVLKTIDEGDADEVTKQRIHDGKEKPGALWHI" +
       "YAAKDAEKIRELLRKVGEEQGQENPPDHDPIHDQSWYLDQTLRKRLYE" +
       "EYGVQGWAIVQFLGDAVFIPAGAPHQVHNLYSCIKVAEDFVSPEHVKH" +
       "CFRLTQEF",
       "-MIPHSWICEKHILWLKDYKNSSNWKLFKECWKQGQPAVVSGVHKKMN" +
       "ISLWKAESISLDFGDHQADLLNCKD-SIISNANVKEFWDGFEEVSKR-" +
       "-----ETVVLKLKDWPSGEDFKTMMPARYEDLLKSLPLPEYCNPEGKF" +
       "NLASHLPGFFVR---GPRLCSAYGVVAAKDHDIGTTNLHIEVSDVVNI" +
       "LVYVGIAKGNGILSKAGILKKFEEEDLDDILRKRLKDSSEIPGALWHI" +
       "YAGKDVDKIREFLQKISKEQG------HDPIRDQSWYVNKKLRQRLLE" +
       "EYGVRTCTLIQFLGDAIVLPAGALHQVQNFHSCIQVTEDFVSPEHLVE" +
       "SFHLTQEL"
    )

    range_ = SequenceRange(48, 96, alignment.target_alignment.replace('-',''))

    template_seq = domain_aligner._get_template_sequence_in_target_range(alignment, range_)
    eq_(template_seq, "ISLWKAESISLDFGDHQADLLNCKD-SIISNANVKEFWDGFEEVSKR-")


def test_find_shared_hits_ranges():
    d = {
        SequenceRange(1, 22, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"):
            TargetTemplateAlignment("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                    "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV"),
        SequenceRange(3, 18, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"):
            TargetTemplateAlignment("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                    "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV"),
        SequenceRange(3, 19, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"):
            TargetTemplateAlignment("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                    "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV")
    }
    template_id = TemplateID('1xxx', 'A')
    for v in d.values():
        v.template_id = template_id

    r = domain_aligner._find_shared_hits_ranges(d)

    eq_(len(r), 1)
    eq_(len(list(r.values())[0]), 3)


def test_remove_enclosing():
    d = {
        SequenceRange(1, 22, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"):
            TargetTemplateAlignment("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                    "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV"),
        SequenceRange(3, 18, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"):
            TargetTemplateAlignment("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                    "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV"),
        SequenceRange(3, 20, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"):
            TargetTemplateAlignment("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                    "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV")
    }
    e = SequenceRange(2, 19, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")

    r = domain_aligner._remove_enclosing(e, d)
    eq_(len(r), 2)


def test_alignment_ok_for_range():
    r = SequenceRange(1, 22, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    a = TargetTemplateAlignment("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                "-----AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA--")

    ok_(domain_aligner._alignment_ok_for_range(r, a))

    r = SequenceRange(1, 22, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    a = TargetTemplateAlignment("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                "-------GGGG------------------------------------")

    ok_(not domain_aligner._alignment_ok_for_range(r, a))


def test_remove_duplicate_ranges():
    r = domain_aligner._remove_duplicate_ranges([
        SequenceRange(1, 22, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
        SequenceRange(1, 22, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    ])
    eq_(len(r), 1)


@with_setup(setup, teardown)
def test_merge_similar_ranges():
    rs = domain_aligner._merge_similar_ranges([
        SequenceRange(1, 22, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
        SequenceRange(19, 30, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
        SequenceRange(0, 23, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
        SequenceRange(1, 23, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    ])

    eq_(len(rs), 2)


def test_is_better_than():
    hit1 = BlastAlignment('', '', '',
                          1, 10, 'AAAAAAAAAA',
                          1, 10, 'AAAAAAAAAA')
    hit2 = BlastAlignment('', '', '',
                          1, 10, 'AAAAAAAAAA',
                          1, 10, '-AAAAAAAAA')
    ok_(domain_aligner._is_better_than(hit1, hit2))
    ok_(not domain_aligner._is_better_than(hit2, hit1))

    hit1 = BlastAlignment('', '', '',
                          1, 10, 'AAAAAAAAAA',
                          1, 10, 'AAAAAAAAAA')
    hit2 = BlastAlignment('', '', '',
                          1, 10, 'AAAAAAAAAA',
                          1, 10, 'VAAAAAAAAA')
    ok_(domain_aligner._is_better_than(hit1, hit2))
    ok_(not domain_aligner._is_better_than(hit2, hit1))


@with_setup(setup, teardown)
def test_filter_forbidden_ranges():
    il = [
        SequenceRange(1, 15, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
        SequenceRange(0, 16, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
        SequenceRange(17, 40, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
        SequenceRange(20, 33, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
    ]
    from hommod import default_settings as settings
    il[1].ac = settings.FORBIDDEN_INTERPRO_DOMAINS[0]

    rs = domain_aligner._filter_forbidden_ranges(il)

    eq_(len(rs), 2)
