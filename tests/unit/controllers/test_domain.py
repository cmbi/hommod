from nose.tools import eq_

from hommod.controllers.domain import domain_aligner
from hommod.models.range import SequenceRange
from hommod.models.align import TargetTemplateAlignment
from hommod.models.template import TemplateID


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
    eq_(len(r.values()[0]), 3)
