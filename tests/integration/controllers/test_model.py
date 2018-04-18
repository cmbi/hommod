import logging

from nose.tools import with_setup, eq_, ok_

from hommod.controllers.blast import blaster
from hommod.services.dssp import dssp
from hommod.controllers.blacklist import blacklister
from hommod.services.interpro import interpro
from hommod.controllers.domain import domain_aligner
from hommod.controllers.model import modeler
from hommod.models.template import TemplateID
from hommod.models.range import SequenceRange
from hommod.models.align import DomainAlignment
from hommod.controllers.clustal import clustal_aligner
from hommod.controllers.kmad import kmad_aligner
import hommod.default_settings as settings


_log = logging.getLogger(__name__)


def setup():
    blaster.blastp_exe = settings.BLASTP_EXE
    dssp.dssp_dir = settings.DSSP_DIR
    blacklister.file_path = settings.BLACKLIST_FILE_PATH
    interpro.url = settings.INTERPRO_URL
    modeler.yasara_dir = settings.YASARA_DIR
    modeler.uniprot_databank = settings.UNIPROT_BLAST_DATABANK
    domain_aligner.forbidden_interpro_domains = settings.FORBIDDEN_INTERPRO_DOMAINS
    domain_aligner.similar_ranges_min_overlap_percentage = settings.SIMILAR_RANGES_MIN_OVERLAP_PERCENTAGE
    domain_aligner.similar_ranges_max_length_difference_percentage = settings.SIMILAR_RANGES_MAX_LENGTH_DIFFERENCE_PERCENTAGE
    domain_aligner.min_percentage_coverage = settings.DOMAIN_MIN_PERCENTAGE_COVERAGE
    domain_aligner.template_blast_databank = settings.TEMPLATE_BLAST_DATABANK
    domain_aligner.max_merge_distance = settings.DOMAIN_MAX_MERGE_DISTANCE
    domain_aligner.highly_homologous_percentage_identity = settings.HIGHLY_HOMOLOGOUS_PERCENTAGE_IDENTITY
    clustal_aligner.clustalw_exe = settings.CLUSTALW_EXE
    kmad_aligner.kmad_exe = settings.KMAD_EXE


def end():
    pass


@with_setup(setup, end)
def test_no_alignment_flip():
    seq = ( 
"MGKLVALVLLGVGLSLVGEMFLAFRERVNASREVEPVEPENCHLIEELESGSEDIDILPSGLAFISSGLKYP" +
"GMPNFAPDEPGKIFLMDLNEQNPRAQALEISGGFDKELFNPHGISIFIDKDNTVYLYVVNHPHMKSTVEIFK" +
"FEEQQRSLVYLKTIKHELLKSVNDIVVLGPEQFYATRDHYFTNSLLSFFEMILDLRWTYVLFYSPREVKVVA" +
"KGFCSANGITVSADQKYVYVADVAAKNIHIMEKHDNWDLTQLKVIQLGTLVDNLTVDPATGDILAGCHPNPM" +
"KLLNYNPEDPPGSEVLRIQNVLSEKPRVSTVYANNGSVLQGTSVASVYHGKILIGTVFHKTLYCEL")

    species_id = 'human'

    range_ = SequenceRange(183, 265, seq)
    template_id = TemplateID('4zrn', 'A')
    alignment = DomainAlignment(
"YFTNSLLSFFEMILDLRWT---YVLFYSPRE-----VKVVA---KGFCSANGITVSAD-Q--K-YVYVADVAAKNIHIMEKHDNWDLTQLKVIQLGT",
"YSTEMYLEFFAREYGLKYTVLRYANVYGPRQDPYGEAGVVAIFTERMLRGEEVHIFGDGEYVRDYVYVDDVVRANLLAMEKGDN------EVFNIGT",
                        range_, template_id)

    context = modeler._prepare_context(alignment.template_id.pdbid)
    context.set_main_target(seq, species_id, alignment.template_id.chain_id)

    chain_alignments = modeler._make_alignments(seq, species_id, alignment, context)
    for chain_id in chain_alignments:

        _log.debug("got alignment {}: {}".format(chain_id, chain_alignments[chain_id]))
        ok_(chain_alignments[chain_id].target_alignment.replace('-','') in seq)


@with_setup(setup, end)
def test_init_template_5GOX():
    context = modeler._prepare_context('5GOX')

    eq_(len(context.get_chain_ids()), 2)


@with_setup(setup, end)
def test_init_template_5MHF():
    context = modeler._prepare_context('5MHF')

    eq_(len(context.get_chain_ids()), 4)
