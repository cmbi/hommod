import os
from hashlib import md5
from glob import glob
import tarfile
from filelock import FileLock

from hommod.models.align import Alignment
from hommod.models.error import InitError
from hommod.controllers.fasta import parse_fasta_from_string


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
        return paths

    def list_models(self, target_sequence, species_id, required_resnum=None, template_id=None):
        if self.model_dir is None:
            raise InitError("model directory is not set")

        sequence_id = self.get_sequence_id(target_sequence)

        species_id = species_id.upper()

        if template_id is None:
            wildcard = "%s_%s_*.tgz" % (sequence_id, species_id)
        else:
            wildcard = "%s_%s_*_%s-%s.tgz" % (sequence_id, species_id, template_id.pdbid, template_id.chain_id)

        wildcard = os.path.join(self.model_dir, wildcard)

        paths = glob(wildcard)

        if required_resnum is None:
            return paths
        else:
            matching_paths = []
            for path in paths:
                name = os.path.splitext(os.path.basename(path))[0]
                range_ = name.split('_')[2]

                start, end = range_.split('-')
                start = int(start)
                end = int(end)
                if required_resnum >= start and required_resnum <= end:
                    matching_paths.append(path)

            return matching_paths

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
                f = ar.extractfile(os.path.join(dir_name, 'align.fa'))
                alignment_fasta = parse_fasta_from_string(f.read())

                rows = {}
                for key in alignment_fasta:
                    rows[key] = alignment_fasta[key].split('|')

                alignments = []
                for n in range(len(rows.values()[0])):
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

                return f.read()
            finally:
                if f is not None:
                    f.close()


model_storage = ModelStorage()
