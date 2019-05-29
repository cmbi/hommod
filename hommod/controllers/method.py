import logging

from hommod.controllers.pdb import parse_chain_order_from_string
from hommod.controllers.clustal import clustal_aligner
from hommod.controllers.storage import model_storage
from hommod.models.align import Alignment


_log = logging.getLogger(__name__)

def select_best_model(tar_paths, sequence, position):

    best_path = None
    best_pid = 0.0
    for tar_path in tar_paths:
        model_name = model_storage.get_model_name_from_path(tar_path)
        sequence_id = model_storage.get_sequence_id_from_name(model_name)

        if position is not None and not model_storage.model_covers(tar_path, sequence, position):
            _log.debug("{} does not cover {}".format(tar_path, position))
            continue

        chain_alignments = model_storage.extract_alignments(tar_path)

        for alignment in chain_alignments:
            if alignment.get_sequence('target') in sequence:

                keys = list(alignment.aligned_sequences.keys())
                if len(keys) != 2:
                    _log.error("alignment for chain {} in {} has too many keys: {}".format(chain_ids[0], tar_path, keys))
                    continue

                pid = alignment.get_percentage_identity(keys[0], keys[1])
                if pid > best_pid:
                    best_pid = pid
                    best_path = tar_path

    return best_path


def select_best_domain_alignment(domain_alignments):
    best = None
    for alignment in domain_alignments:
        if best is None or alignment.get_percentage_identity() > best.get_percentage_identity():
            best = alignment

    return best
