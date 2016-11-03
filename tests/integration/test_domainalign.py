from hommod_rest.services.domainalign import getAlignments
from nose.tools import ok_, eq_, with_setup

from hommod_rest.services.interpro import interpro
from hommod_rest.services.blast import blaster
from hommod_rest.services.secstr import secstr
from hommod_rest.services.align import aligner
from hommod_rest.services.model import modeler

import hommod_rest.default_settings as config

import os

mydir = os.path.dirname(__file__)


def setupApps():
    interpro.interpro_exe = config.INTERPROSCAN
    interpro.storage_dir = config.INTERPRODIR
    blaster.uniprot_db = os.path.join(mydir, "data/mini")  # small version
    blaster.templates_db = config.TEMPLATESDB
    blaster.blastp_exe = config.BLASTP
    secstr.dssp_dir = config.DSSPDIR
    secstr.yasara_dir = config.YASARADIR
    aligner.clustal_exe = config.CLUSTAL
    aligner.kmad_exe = config.KMAD
    modeler.yasara_dir = config.YASARADIR
    modeler.execution_root_dir = config.EXECUTIONDIR
    modeler.model_root_dir = config.MODELDIR
    modeler.template_blacklist = config.TEMPLATE_BLACKLIST


def tearDown():
    pass


@with_setup(setupApps, tearDown)
def test_1CBO_A():
    seq = "GYVPAVVIGTGYGAAVSALRLGEAGVQTLMLEMGQLWNQPGPDGNIFCGMLNPDKRSS" \
        + "WFKNRTEAPLGSFLWLDVVNRNIDPYAGVLDRVNYDQMSVYVGRGVGGGSLVNGGMAV" \
        + "EPKRSYFEEILPRVDSSEMYDRYFPRANSMLRVNHIDTKWFEDTEWYKFARVSREQAG" \
        + "KAGLGTVFVPNVYDFGYMQREAAGEVPKSALATEVIYGNNHGKQSLDKTYLAAALGTG" \
        + "KVTIQTLHQVKTIRQTKDGGYALTVEQKDTDGKLLATKEISCRYLFLGAGSLGSTELL" \
        + "VRARDTGTLPNLNSEVGAGWGPNGNIMTARANHMWNPTGAHQSSIPALGIDAWDNSDS" \
        + "SVFAEIAPMPAGLETWVSLYLAITKNPQRGTFVYDAATDRAKLNWTRDQNAPAVNAAK" \
        + "ALFDRINKANGTIYRYDLFGTQLKAFADDFCYNPLGGCVLGKATDDYGRVAGYKNLYV" \
        + "TDGSLIPGSVGVNPFVTITALAERNVERIIKQDV"

    domainranges = interpro.get_domain_locations(seq)
    ok_(len(domainranges) > 0)

    alignments = getAlignments(domainranges, seq)

    eq_(type(alignments), list)
    ok_(len(alignments) > 0)
