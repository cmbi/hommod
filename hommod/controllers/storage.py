import os
import logging
from hashlib import md5
from glob import glob
import tarfile
from filelock import FileLock

from hommod.models.align import Alignment
from hommod.models.error import InitError
from hommod.controllers.fasta import parse_fasta_from_string
from hommod.controllers.pdb import parse_seqres_from_string
from hommod.controllers.clustal import clustal_aligner


_log = logging.getLogger(__name__)


class ModelStorage:
    def __init__(self, model_dir=None):
        self.model_dir = model_dir

    def get_sequence_id(self, sequence):
        hash_ = md5(sequence.encode('ascii')).hexdigest()
        return hash_

    def list_all_models(self):
        if self.model_dir is None:
            raise InitError("model directory is not set")

        wildcard = os.path.join(self.model_dir, "*.tgz")

        paths = glob(wildcard)
        paths = [path for path in paths if '_error' not in path]
        return paths

    def list_models(self, target_sequence, species_id, required_resnum=None, template_id=None):
        if self.model_dir is None:
            raise InitError("model directory is not set")
        elif not os.path.isdir(self.model_dir):
            raise InitError("No such directory: {}".format(self.model_dir))


        sequence_id = self.get_sequence_id(target_sequence)

        species_id = species_id.upper()

        if template_id is None:
            wildcard = "%s_%s_*.tgz" % (sequence_id, species_id)
        else:
            case_insensitive_pdbid = ""
            for i in range(len(template_id.pdbid)):
                char = template_id.pdbid[i]
                if char.isalpha():
                    case_insensitive_pdbid += '[%s%s]' % (char.lower(), char.upper())
                else:
                    case_insensitive_pdbid += char

            wildcard = "%s_%s_*_%s-%s.tgz" % (sequence_id, species_id, case_insensitive_pdbid, template_id.chain_id)

        wildcard = os.path.join(self.model_dir, wildcard)

        paths = glob(wildcard)
        paths = [path for path in paths if '_error' not in path]

        if required_resnum is None:
            return paths
        else:
            matching_paths = []
            for path in paths:
                if self.model_covers(path, target_sequence, required_resnum):
                    matching_paths.append(path)

            return matching_paths

    def model_covers(self, tar_path, sequence, covered_residue_number):
        sequence_id = self.get_sequence_id(sequence)

        seqres_sequences = parse_seqres_from_string(self.extract_model(tar_path))
        _log.debug(str(seqres_sequences))

        for chain_id in seqres_sequences:
            _log.debug(chain_id)

            if len(seqres_sequences[chain_id]) <= 0:
                continue

            full_to_model = clustal_aligner.align({'model': ''.join([aa.letter for aa in seqres_sequences[chain_id]]),
                                                   'full': sequence})
            _log.debug(str(full_to_model))

            resnum = 0
            for i in range(len(full_to_model.aligned_sequences['full'])):
                if full_to_model.aligned_sequences['full'][i].isalpha():
                    resnum += 1

                    if resnum == covered_residue_number and full_to_model.aligned_sequences['model'][i].isalpha():
                        return True
        return False

    def get_model_name(self, main_target_sequence, target_species_id,
                       main_domain_alignment, template_id):

        if template_id is None:
            name = "%s_%s_%i-%i" % (self.get_sequence_id(main_target_sequence),
                                    target_species_id.upper(),
                                    main_domain_alignment.range.start + 1, main_domain_alignment.range.end)
        else:
            name = "%s_%s_%i-%i_%s" % (self.get_sequence_id(main_target_sequence),
                                       target_species_id.upper(),
                                       main_domain_alignment.range.start + 1, main_domain_alignment.range.end,
                                       str(template_id))
        return name

    def get_sequence_id_from_name(self, model_name):
        s = model_name.split('_')
        return s[0]

    def get_model_name_from_path(self, tar_path):
        return os.path.splitext(os.path.basename(tar_path))[0]

    def get_tar_path_from_name(self, name):
        if self.model_dir is None:
            raise InitError("model directory is not set")

        return os.path.join(self.model_dir, name + '.tgz')

    def get_tar_path(self, target_sequence, target_species_id, main_domain_alignment, template_id):
        name = self.get_model_name(target_sequence, target_species_id,
                                   main_domain_alignment, template_id)

        return self.get_tar_path_from_name(name)

    def get_error_tar_path_from_name(self, name):
        if self.model_dir is None:
            raise InitError("model directory is not set")

        return os.path.join(self.model_dir, name + '_error.tgz')

    def get_error_tar_path(self, target_sequence, target_species_id, main_domain_alignment, template_id):
        name = self.get_model_name(target_sequence, target_species_id,
                                   main_domain_alignment, template_id)

        return self.get_error_tar_path_from_name(name)

    def get_model_lock(self, main_target_sequence, target_species_id,
                       main_domain_alignment, template_id):
        if self.model_dir is None:
            raise InitError("model directory is not set")

        lock_name = 'lock_model_' + self.get_model_name(main_target_sequence,
                                                        target_species_id,
                                                        main_domain_alignment,
                                                        template_id)
        lock_path = os.path.join(self.model_dir, lock_name)
        return FileLock(lock_path)

    def extract_alignments(self, tar_path):
        dir_name = os.path.splitext(os.path.basename(tar_path))[0]
        with tarfile.open(tar_path, 'r:gz') as ar:
            f = None
            try:
                path = os.path.join(dir_name, 'align.fa')
                if path not in ar.getnames():
                    path = os.path.join(dir_name, 'align.fasta')

                f = ar.extractfile(path)
                alignment_fasta = parse_fasta_from_string(f.read().decode('ascii'))

                rows = {}
                for key in alignment_fasta:
                    rows[key] = alignment_fasta[key].split('|')

                alignments = []
                for n in range(len(list(rows.values())[0])):
                    a = {key: rows[key][n] for key in rows}
                    alignments.append(Alignment(a))
                return alignments
            finally:
                if f is not None:
                    f.close()

    def extract_selected_targets(self, tar_path):
        dir_name = os.path.splitext(os.path.basename(tar_path))[0]
        with tarfile.open(tar_path, 'r:gz') as ar:
            f = None
            try:
                targets = {}

                f = ar.extractfile(os.path.join(dir_name, 'selected-targets.txt'))
                for line in f:
                    line = line.decode('ascii')
                    if ':' in line:
                        chain_id, target_id = line.split(':')
                        targets[chain_id.strip()] = target_id.strip()

                return targets
            finally:
                if f is not None:
                    f.close()

    def extract_model(self, tar_path):
        dir_name = os.path.splitext(os.path.basename(tar_path))[0]
        with tarfile.open(tar_path, 'r:gz') as ar:
            f = None
            try:
                f = ar.extractfile(os.path.join(dir_name, 'target.pdb'))

                return f.read().decode('ascii')
            finally:
                if f is not None:
                    f.close()


model_storage = ModelStorage()
