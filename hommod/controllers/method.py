import logging

from hommod.controllers.pdb import parse_chain_order_from_string
from hommod.controllers.clustal import clustal_aligner
from hommod.controllers.storage import model_storage
from hommod.models.align import Alignment


_log = logging.getLogger(__name__)

def select_best_model(tar_paths):

    best_path = None
    best_pid = 0.0
    for tar_path in tar_paths:
        model_name = model_storage.get_model_name_from_path(tar_path)
        sequence_id = model_storage.get_sequence_id_from_name(model_name)

        pdb_str = model_storage.extract_model(tar_path)
        chain_order = parse_chain_order_from_string(pdb_str)

        chain_alignments = model_storage.extract_alignments(tar_path)

        chain_targets = model_storage.extract_selected_targets(tar_path)
        chain_ids = list(filter(lambda chain_id: chain_targets[chain_id] == sequence_id,
                                chain_targets.keys()))
        if len(chain_ids) <= 0:
            _log.error("no main target {} in selected targets of {}".format(sequence_id, tar_path))
            continue

        chain_i = chain_order.index(chain_ids[0])
        if chain_i == -1:
            _log.error("chain {} not in {}".format(chain_ids[0], tar_path))
            continue

        alignment = chain_alignments[chain_i]
        keys = alignment.aligned_sequences.keys()
        if len(keys) != 2:
            _log.error("alignment for chain {} in {} has too many keys: {}".format(chain_ids[0], tar_path, keys))
            continue

        pid = alignment.get_percentage_identity(keys[0], keys[1])
        if pid > best_pid:
            best_pid = pid
            best_path = tar_path

    return best_path
