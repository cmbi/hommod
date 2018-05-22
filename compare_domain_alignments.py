import random
from requests.exceptions import HTTPError

from hommod.default_settings import (SPROT_FASTA, DSSP_DIR, INTERPRO_URL,
                                     HIGHLY_HOMOLOGOUS_PERCENTAGE_IDENTITY,
                                     DOMAIN_MAX_MERGE_DISTANCE,
                                     DOMAIN_MIN_PERCENTAGE_COVERAGE,
                                     SIMILAR_RANGES_MIN_OVERLAP_PERCENTAGE,
                                     SIMILAR_RANGES_MAX_LENGTH_DIFFERENCE_PERCENTAGE,
                                     FORBIDDEN_INTERPRO_DOMAINS,
                                     TEMPLATE_BLAST_DATABANK, KMAD_EXE, BLASTP_EXE,
                                     CACHE_REDIS_HOST, CACHE_REDIS_PORT, CACHE_REDIS_DB,
                                     CACHE_EXPIRATION_TIME, CACHE_LOCK_TIMEOUT,
                                     BLACKLIST_FILE_PATH)
from hommod.controllers.fasta import parse_fasta
from hommod.controllers.domain import domain_aligner
from hommod.controllers.blast import blaster
from hommod.controllers.kmad import kmad_aligner
from hommod.services.dssp import dssp
from hommod.services.interpro import interpro
from hommod.services.helpers.cache import cache_manager as cm
from hommod.controllers.blacklist import blacklister


blacklister.file_path = BLACKLIST_FILE_PATH
cm.redis_hostname = CACHE_REDIS_HOST
cm.redis_port = CACHE_REDIS_PORT
cm.redis_db = CACHE_REDIS_DB
cm.expiration_time = CACHE_EXPIRATION_TIME
cm.lock_timeout = CACHE_LOCK_TIMEOUT
dssp.dssp_dir = DSSP_DIR
interpro.url = INTERPRO_URL
domain_aligner.highly_homologous_percentage_identity = HIGHLY_HOMOLOGOUS_PERCENTAGE_IDENTITY
domain_aligner.max_merge_distance = DOMAIN_MAX_MERGE_DISTANCE
domain_aligner.min_percentage_coverage = DOMAIN_MIN_PERCENTAGE_COVERAGE
domain_aligner.forbidden_interpro_domains = FORBIDDEN_INTERPRO_DOMAINS
domain_aligner.template_blast_databank = TEMPLATE_BLAST_DATABANK
domain_aligner.similar_ranges_min_overlap_percentage = SIMILAR_RANGES_MIN_OVERLAP_PERCENTAGE
domain_aligner.similar_ranges_max_length_difference_percentage = SIMILAR_RANGES_MAX_LENGTH_DIFFERENCE_PERCENTAGE
kmad_aligner.kmad_exe = KMAD_EXE
blaster.blastp_exe = BLASTP_EXE


def pick_random_sequences(n):
    sprot_sequences = parse_fasta(SPROT_FASTA)
    keys = random.sample(sprot_sequences.keys(), n)

    return {key:sprot_sequences[key] for key in keys}


sequences = pick_random_sequences(10)
for key in sequences:

    while True:
        try:
            domain_alignments = domain_aligner.get_domain_alignments(sequences[key])
            break
        except HTTPError:
            continue

    for domain_alignment in domain_alignments:
        template_seq = dssp.get_sequence(domain_alignment.template_id)
        template_secstr = dssp.get_secondary_structure(domain_alignment.template_id)
        full_alignment = kmad_aligner.align(template_seq, template_secstr, sequences[key])

        print(key, domain_alignment.template_id, domain_alignment.get_percentage_identity(),
                                                 full_alignment.get_percentage_identity())
