import tempfile
import os
import sys
import shutil
from multiprocessing import Process, Value
from argparse import ArgumentParser
import logging

from hommod.controllers.fasta import parse_fasta
from hommod.models.template import TemplateID
from hommod.controllers.domain import domain_aligner
from hommod.controllers.model import modeler
from hommod.controllers.soup import soup
from hommod.services.uniprot import uniprot
from hommod.services.interpro import interpro
from hommod.controllers.kmad import kmad_aligner
from hommod.controllers.clustal import clustal_aligner
from hommod.controllers.storage import model_storage
from hommod.controllers.blast import blaster
from hommod.controllers.blacklist import blacklister
from hommod.services.dssp import dssp
from hommod.services.helpers.cache import cache_manager as cm
import hommod.default_settings as settings


class ModelProcess(Process):
    def __init__(self, sequence, species_id, domain_alignment, output_dir):
        Process.__init__(self)
        self.daemon = True

        self.exc_info = None

        self.sequence = sequence
        self.species_id = species_id
        self.domain_alignment = domain_alignment
        self.output_dir = output_dir

    def run(self):
        path = modeler.build_model(self.sequence, self.species_id, self.domain_alignment)

        _log.info(path)
        shutil.copy(path, self.output_dir)

        soup.yasara.Exit()


logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
_log = logging.getLogger(__name__)

if __name__ == "__main__":

    cm.disable()
    uniprot.url = settings.UNIPROT_URL
    interpro.url = settings.INTERPRO_URL
    kmad_aligner.kmad_exe = settings.KMAD_EXE
    clustal_aligner.clustalw_exe = settings.CLUSTALW_EXE

    soup.yasara_dir = settings.YASARA_DIR

    modeler.uniprot_databank = settings.UNIPROT_BLAST_DATABANK

    domain_aligner.forbidden_interpro_domains = settings.FORBIDDEN_INTERPRO_DOMAINS
    domain_aligner.similar_ranges_min_overlap_percentage = settings.SIMILAR_RANGES_MIN_OVERLAP_PERCENTAGE
    domain_aligner.similar_ranges_max_length_difference_percentage = settings.SIMILAR_RANGES_MAX_LENGTH_DIFFERENCE_PERCENTAGE
    domain_aligner.min_percentage_coverage = settings.DOMAIN_MIN_PERCENTAGE_COVERAGE
    domain_aligner.template_blast_databank = settings.TEMPLATE_BLAST_DATABANK
    domain_aligner.max_merge_distance = settings.DOMAIN_MAX_MERGE_DISTANCE
    domain_aligner.highly_homologous_percentage_identity = settings.HIGHLY_HOMOLOGOUS_PERCENTAGE_IDENTITY

    blaster.blastp_exe = settings.BLASTP_EXE
    blacklister.file_path = settings.BLACKLIST_FILE_PATH
    dssp.dssp_dir = settings.DSSP_DIR


    arg_parser = ArgumentParser(description="Build Models for given criteria")
    arg_parser.add_argument('--output-dir', help="output dir to put the models in")
    arg_parser.add_argument('fasta', help="fasta with input target sequence")
    arg_parser.add_argument('species', help="target species id")
    arg_parser.add_argument('--position', help="residue position that the models should cover", type=int)
    arg_parser.add_argument('--template', help="underscore separated template pdbid and chain")

    args = arg_parser.parse_args()


    tmp_dir = tempfile.mkdtemp()
    model_storage.model_dir = tmp_dir

    final_output_dir = settings.MODEL_DIR
    if args.output_dir:
        final_output_dir = args.output_dir

    if not os.path.isdir(final_output_dir):
        raise ValueError("Not a directory: {}".format(final_output_dir))

    try:
        sequence = parse_fasta(args.fasta).values()[0]

        species_id = args.species.upper()

        if args.template:
            pdbid, chain_id = args.template.split('_')
            template_id = TemplateID(pdbid, chain_id)
        else:
            template_id = None

        domain_alignments = domain_aligner.get_domain_alignments(sequence, args.position, template_id)
        _log.info("{} domain alignments".format(len(domain_alignments)))

        ps = [ModelProcess(sequence, species_id, ali, final_output_dir) for ali in domain_alignments]
        for p in ps:
            p.start()
        for p in ps:
            p.join()
    finally:
        shutil.rmtree(tmp_dir)
