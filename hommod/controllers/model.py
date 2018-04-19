import shutil
import os
import tarfile
import logging
import tempfile
import re

from hommod.controllers.blast import blaster
from hommod.controllers.fasta import write_fasta, parse_fasta_from_string
from hommod.controllers.kmad import kmad_aligner
from hommod.controllers.clustal import clustal_aligner
from hommod.controllers.domain import domain_aligner
from hommod.controllers.blacklist import blacklister
from hommod.models.align import TargetTemplateAlignment, Alignment
from hommod.models.template import TemplateID
from hommod.controllers.soup import soup
from hommod.models.error import TemplateError, ModelRunError, InitError
from hommod.services.pdb import get_pdb_contents
from hommod.services.uniprot import uniprot
from hommod.controllers.storage import model_storage
from hommod.controllers.sequence import is_amino_acid_char


_log = logging.getLogger(__name__)


def swap(x, y):
    return (y, x)


class Modeler:

    def __init__(self, uniprot_databank=None):

        self.uniprot_databank = uniprot_databank

    def build_model(self, main_target_sequence, target_species_id, main_domain_alignment):

        tar_path = model_storage.get_tar_path(main_target_sequence,
                                              target_species_id,
                                              main_domain_alignment,
                                              main_domain_alignment.template_id)

        with model_storage.get_model_lock(main_target_sequence, target_species_id,
                                          main_domain_alignment, main_domain_alignment.template_id):
            if not os.path.isfile(tar_path):

                self._prepare_soup(main_domain_alignment.template_id.pdbid)

                try:
                    # If the template is the same as the target, do no modeling:
                    if main_domain_alignment.get_template_sequence() == soup.get_sequence(main_domain_alignment.template_id.chain_id) and \
                            main_domain_alignment.get_percentage_identity() >= 100.0:

                        main_domain_alignment.target_id = model_storage.get_sequence_id(main_target_sequence)

                        tar_path = self._wrap_template(main_target_sequence, target_species_id,
                                                       main_domain_alignment, main_domain_alignment.template_id)
                        return tar_path


                    soup.set_main_target(main_target_sequence, target_species_id,
                                         main_domain_alignment.template_id.chain_id)

                    chain_alignments = self._make_alignments(main_target_sequence, target_species_id,
                                                             main_domain_alignment)

                    # Delete chains that aren't in the alignment set:
                    for chain_id in soup.get_chain_ids():
                        if chain_id not in chain_alignments:
                            self._delete_chain(chain_id)

                    _log.debug("final alignments: {}".format([(chain_id, chain_alignments[chain_id])
                                                              for chain_id in soup.get_chain_ids()]))
                    _log.debug("final template {} {}".format(soup.template_pdbid,
                                                             [(chain_id, soup.get_sequence(chain_id))
                                                              for chain_id in soup.get_chain_ids()]))

                    tar_path = self._model_run(main_domain_alignment, chain_alignments)
                finally:
                    soup.cleanup()

        return tar_path

    def _prepare_soup(self, template_pdbid):
        self._init_template(template_pdbid)
        try:
            self._oligomerize_template()
        except:
            self._init_template(template_pdbid)

        try:
            self._build_template_symmetry_residues()
        except:
            pass

        self._delete_solvent_residues()
        self._fix_template_errors()

        soup.yasara.CleanObj(soup.template_obj)

    def _init_template(self, template_pdbid):
        soup.template_pdbid = template_pdbid
        soup.yasara.Clear()
        soup.template_obj = soup.yasara.LoadPDB(template_pdbid, download='yes')[0]
        soup.yasara.DelObj("not %i" % soup.template_obj)

        chain_ids = soup.get_chain_ids()
        if len(chain_ids) <= 0:
            raise ValueError("No protein chains in %s" % template_pdbid)

    def _oligomerize_template(self):
        count_chains_unoligomerized = len(soup.get_chain_ids())

        oligo_list = soup.yasara.OligomerizeObj(soup.template_obj)
        if len(oligo_list) > 1:
            for obj in oligo_list:
                if obj != soup.template_obj:
                    soup.yasara.JoinObj(obj, soup.template_obj)

        count_chains_oligomerized = len(soup.get_chain_ids())

        if count_chains_oligomerized < count_chains_unoligomerized:
            raise TemplateError("oligomerisation of %s reduced the number of chains" % soup.template_pdbid)

    def _build_template_symmetry_residues(self):
        soup.yasara.BuildSymRes('obj %i' % soup.template_obj)

    def _delete_solvent_residues(self):
        soup.yasara.DelRes('HOH H2O DOD D2O TIP WAT SOL ACT ACM ACY ' +
                           'EDO EOH FLC FMT GOL I3M IPA MLI MOH PEO ' +
                           'PO4 SO3 SO4 _TE UNX ACE')

    def _fix_template_errors(self):
        self._join_duplicant_chains()
        self._swap_problem_residues()

    def _join_duplicant_chains(self):
        chain_ids = soup.get_chain_ids()
        chain_id_count = {}
        for chain_id in chain_ids:
            if chain_id not in chain_id_count:
                chain_id_count[chain_id] = 0
            chain_id_count[chain_id] += 1

        for chain_id in chain_id_count:
            if chain_id_count[chain_id] > 1:
                rs = soup.yasara.ListMol("obj %i mol %s and protein" % (soup.template_obj, chain_id))
                _log.debug("chain {} listmol: {}".format(chain_id, rs))
                for r in rs[1:]:
                    atom_nr = int(r.split()[-1])
                    soup.yasara.JoinMol("obj %i atom %i" % (soup.template_obj, atom_nr))

                    # Add the peptide bond:
                    soup.yasara.AddBond("obj %i and res atom %i and atom C" % (soup.template_obj, atom_nr - 1),
                                        "obj %i and res atom %i and atom N" % (soup.template_obj, atom_nr))

    def _swap_problem_residues(self):
        soup.yasara.SwapRes('Protein and UNK and atom CA and atom CB', 'ALA')
        soup.yasara.SwapRes('Protein and UNK and atom CA and not atom CB', 'GLY')
        soup.yasara.SwapRes('Protein and CAS', 'CYS')

    def _group_identical_chains(self):

        sequences = {chain_id: soup.get_sequence(chain_id)
                     for chain_id in soup.get_chain_ids()}

        # If there's only 1 chain, then we don't need to do anything:
        if len(sequences) <= 1:
            return sequences.keys()

        alignments = {} 

        grouped = [] 

        ids = sequences.keys()
        while len(ids) > 0: 
            id_ = ids[0]
            ids.remove(id_)
            grouped.append([id_])

            for other_id in ids[:]:

                pair = (id_, other_id)

                # Aligning them all in one run can be very time-consuming,
                # so align two at the time:
                if pair not in alignments:
                    alignments[pair] = clustal_aligner.align({id_: sequences[id_],
                                                              other_id: sequences[other_id]})

                if alignments[pair].get_percentage_identity(id_, other_id) >= 99.0:
                    grouped[-1].append(other_id)
                    ids.remove(other_id)

        return grouped

    def _pick_template_chains(self, main_chain_id, main_target_sequence):

        identical_chain_groups = self._group_identical_chains()
        for group in identical_chain_groups:
            if main_chain_id in group:
                return group

        raise ModelRunError("chain not found in identical groups: {}".format(main_chain_id))

    def _make_alignments(self, main_target_sequence, target_species_id,
                         main_domain_alignment):
        alignments = {}

        # Choose what chains to align the main_target_on
        main_target_chain_ids = self._pick_template_chains(main_domain_alignment.template_id.chain_id,
                                                           main_target_sequence)

        for chain_id in main_target_chain_ids:

            template_chain_sequence = soup.get_sequence(chain_id)
            template_chain_secstr = soup.get_secondary_structure(chain_id)

            alignments[chain_id] = kmad_aligner.align(template_chain_sequence, template_chain_secstr,
                                                      main_domain_alignment.get_target_sequence())
            alignments[chain_id].target_id = model_storage.get_sequence_id(main_target_sequence)

        # Try to find and align target sequences for interacting chains in the template,
        # while keeping in mind which residues interact and must thus be covered by the alignment.
        # We expand the set of involved template chains with every iteration,
        # until all template chains have been added.
        while len(alignments) < len(soup.get_chain_ids()):

            # First, make python remember to which chains the candidate chains interact:
            candidate_chains_interacts_with = {}
            for aligned_chain_id in alignments:
                for interacting_chain_id in soup.list_interacting_chains(aligned_chain_id):
                    # Skip those that we've already aligned, to prevent infinite loops:
                    if interacting_chain_id in alignments:
                        continue

                    if interacting_chain_id not in candidate_chains_interacts_with:
                        candidate_chains_interacts_with[interacting_chain_id] = []
                    candidate_chains_interacts_with[interacting_chain_id].append(aligned_chain_id)

            if len(candidate_chains_interacts_with) <= 0:
                break  # Nothing more to add

            # iterate over chains that might interact with the chains that are already in the set:
            for candidate_chain_id in candidate_chains_interacts_with:

                interacting_chain_alignments = {interacting_chain_id: alignments[interacting_chain_id]
                                                for interacting_chain_id in candidate_chains_interacts_with[candidate_chain_id]}

                template_chain_sequence = soup.get_sequence(candidate_chain_id)
                template_chain_secstr = soup.get_secondary_structure(candidate_chain_id)

                potential_target_sequences = self._find_target_sequences(template_chain_sequence,
                                                                         target_species_id)

                alignments[candidate_chain_id] = self._choose_best_target_alignment(interacting_chain_alignments,
                                                                                    potential_target_sequences,
                                                                                    candidate_chain_id)
                if alignments[candidate_chain_id] is None:
                    alignments[candidate_chain_id] = self._make_poly_A_alignment(candidate_chain_id)
                    alignments[candidate_chain_id].target_id = "poly-A"

        return alignments

    def _find_target_sequences(self, template_chain_sequence, target_species_id):
        if self.uniprot_databank is None:
            raise InitError("species databank dir not set")

        target_sequences = {}

        hits = blaster.blastp(template_chain_sequence, self.uniprot_databank)
        for hit_id in hits:
            if not hit_id.endswith('_' + target_species_id.upper()):
                continue

            for alignment in hits[hit_id]:
                ac = alignment.get_hit_accession_code()
                pid = alignment.get_percentage_identity()
                pcov = alignment.get_percentage_coverage()
                if pid > 70.0:
                    if pcov > 90.0:
                        target_sequences[ac] = uniprot.get_sequence(ac)
        return target_sequences

    def _preserves_interactions(self,
                                candidate_target_segment, candidate_chain_id,
                                interacting_chain_alignments):

        # The pdb file in the soup can be different from the blast hit
        # So make an alignment first!
        template_chain_sequence = soup.get_sequence(candidate_chain_id)
        template_chain_secstr = soup.get_secondary_structure(candidate_chain_id)
        candidate_alignment = kmad_aligner.align(template_chain_sequence,
                                                 template_chain_secstr,
                                                 candidate_target_segment)

        candidate_residue_indices = candidate_alignment.get_covered_template_residues_indices()
        candidate_residues = soup.get_residues(candidate_chain_id)
        covered_candidate_residues = [candidate_residues[i]
                                      for i in candidate_residue_indices]

        for chain_id in interacting_chain_alignments:

            covered_template_residue_indices = \
                interacting_chain_alignments[chain_id].get_covered_template_residues_indices()

            chain_residues = soup.get_residues(chain_id)
            covered_residues = [chain_residues[i]
                                for i in covered_template_residue_indices]

            _log.debug("checking chain {} {} residues against chain {} {} residues for interaction"
                       .format(candidate_chain_id, len(covered_candidate_residues),
                               chain_id, len(covered_residues)))

            # Check every target-covered residue.
            # Return True if a single interacting residue pair is found:
            for candidate_residue in covered_candidate_residues:
                if soup.residue_interacts_with(candidate_residue, covered_residues):
                    return True

        return False

    def _make_poly_A_alignment(self, chain_id):

        template_seq = soup.get_sequence(chain_id)

        return TargetTemplateAlignment('A' * len(template_seq), template_seq)

    def _choose_best_target_alignment(self, interacting_chain_alignments,
                                      potential_target_sequences,
                                      chain_id):

        best_alignment = None

        for target_id in potential_target_sequences:
            template_chain_sequence = soup.get_sequence(chain_id)
            template_chain_secstr = soup.get_secondary_structure(chain_id)

            alignment = kmad_aligner.align(template_chain_sequence,
                                           template_chain_secstr,
                                           potential_target_sequences[target_id])

            _log.debug("alignment {} has coverage {} %".format(alignment, alignment.get_percentage_coverage()))

            if alignment.get_percentage_coverage() < 90.0:
                # If the coverage is too low, we need to bother interpro.
                domain_alignments = \
                    domain_aligner.get_domain_alignments(potential_target_sequences[target_id],
                                                         None,
                                                         TemplateID(soup.template_pdbid, chain_id))

                interacting_alignments = list(filter(lambda ali: self._preserves_interactions(ali.target_alignment.replace('-', ''),
                                                                                              chain_id, interacting_chain_alignments),
                                                     domain_alignments))

                _log.debug("preserve interactions with chains {}: filtered {} alignments out of {}"
                           .format(interacting_chain_alignments.keys(),
                                   len(interacting_alignments), len(domain_alignments)))

                if len(interacting_alignments) > 0:
                    domain_alignment = self._join_alignments_to_best_template_coverage(interacting_alignments)
                elif len(domain_alignments) > 0:
                    domain_alignment = self._join_alignments_to_best_template_coverage(domain_alignments)
                else:
                    continue

                alignment = kmad_aligner.align(template_chain_sequence,
                                               template_chain_secstr,
                                               domain_alignment.get_target_sequence())

            alignment.target_id = target_id

            if best_alignment is None or \
                    best_alignment.get_percentage_identity() < alignment.get_percentage_identity():
                best_alignment = alignment

        return best_alignment

    def _join_alignments_to_best_template_coverage(self, domain_alignments):

        # First determine which alignment has the largest cover:
        best_alignment = None
        for alignment in domain_alignments:
            if alignment.template_id != domain_alignments[0].template_id:
                raise ValueError("not all the same template")

            if best_alignment is None or \
                    alignment.get_percentage_coverage() > best_alignment.get_percentage_coverage():
                best_alignment = alignment

        # Now merge the best alignment with smaller, compatible alignments
        for alignment in domain_alignments:
            if self._alignments_compatible(best_alignment, alignment):
                best_alignment = self._merge_alignments(best_alignment, alignment)

        return best_alignment

    def _alignments_compatible(self, alignment1, alignment2):
        span1 = alignment1.get_relative_span()
        span2 = alignment2.get_relative_span()

        # They must not occupy the same piece of template!
        return not span1.overlaps_with(span2)

    def _merge_alignments(self, alignment1, alignment2):
        span1 = alignment1.get_relative_span()
        span2 = alignment2.get_relative_span()

        # make sure span1 is the left sided:
        if span2 < span1:
            span1, span2 = swap(span1, span2)
            alignment1, alignment2 = swap(alignment1, alignment2)

        # Determine alignment positions of span1 and span2:
        i1 = 0
        naa = 0
        while naa < span1.end:
            if is_amino_acid_char(alignment1.template_alignment[i1]):
                naa += 1
            i1 += 1

        i2 = 0
        naa = 0
        while naa < span2.start:
            if is_amino_acid_char(alignment2.template_alignment[i2]):
                naa += 1
            i2 += 1

        return TargetTemplateAlignment(alignment1.target_alignment[: i1] + alignment2.target_alignment[: i2],
                                       alignment1.template_alignment[: i1] + alignment2.template_alignment[: i2])

    def _write_model_alignment_fasta(self, chain_alignments, fasta_path):
        # Creates the input file for YASAYA's modeling run.
        # Includes all the chain ids in 'chain_order' in the alignment.

        chain_order = soup.get_chain_ids()

        with open(fasta_path, 'w') as f:
            f.write('>target\n')
            n = 0
            for chain_id in chain_order:
                if chain_id not in chain_alignments:
                    continue
                if n > 0:
                    f.write('|')
                f.write(chain_alignments[chain_id].target_alignment)
                n += 1
            f.write('\n')

            f.write('>%s\n' % soup.template_pdbid)
            n = 0
            for chain_id in chain_order:
                if chain_id not in chain_alignments:
                    continue
                if n > 0:
                    f.write('|')
                f.write(chain_alignments[chain_id].template_alignment)
                n += 1
            f.write('\n')

    def _wrap_template(self, main_target_sequence, target_species_id, main_domain_alignment, template_id):
        model_name = model_storage.get_model_name(main_target_sequence, target_species_id,
                                                  main_domain_alignment, template_id)

        work_dir_path = tempfile.mkdtemp()
        align_fasta_path = os.path.join(work_dir_path, 'align.fa')

        try:
            os.chdir(work_dir_path)

            write_fasta(align_fasta_path, {'target': main_domain_alignment.target_alignment,
                                           str(template_id): main_domain_alignment.template_alignment})

            model_path = os.path.join(work_dir_path, 'target.pdb')
            with open(model_path, 'w') as f:
                f.write(get_pdb_contents(template_id.pdbid))

            self._write_selected_targets({template_id.chain_id: main_domain_alignment},
                                         os.path.join(work_dir_path, 'selected-targets.txt'))

            tar_path = model_storage.get_tar_path(main_target_sequence,
                                                  target_species_id,
                                                  main_domain_alignment,
                                                  template_id)
            with tarfile.open(tar_path, mode="w:gz") as ar:
                ar.add(work_dir_path, arcname=model_name)

            return tar_path
        finally:
            if os.path.isdir(work_dir_path):
                shutil.rmtree(work_dir_path)

    def _model_run(self, main_domain_alignment, chain_alignments):

        model_name =  model_storage.get_model_name(soup.get_main_target_sequence(),
                                                   soup.target_species_id,
                                                   main_domain_alignment,
                                                   TemplateID(soup.template_pdbid,
                                                              soup.main_target_chain_id))

        work_dir_path = tempfile.mkdtemp()
        align_fasta_path = os.path.join(work_dir_path, 'align.fa')
        output_yob_path = os.path.join(work_dir_path, 'target.yob')
        error_path = os.path.join(work_dir_path, 'errorexit.txt')
        error_scene_path = os.path.join(work_dir_path, 'errorexit.sce')
        try:
            soup.yasara.CD(work_dir_path)

            self._write_model_alignment_fasta(chain_alignments, align_fasta_path)

            soup.yasara.Processors(1)

            soup.yasara.ExperimentHomologyModeling(templateobj=soup.template_obj,
                                                   alignfile=align_fasta_path,
                                                   templates="1, sameseq = 1",
                                                   alignments=1,
                                                   termextension=0,
                                                   oligostate=32,
                                                   looplenmax=10,
                                                   animation='fast',
                                                   speed='fast',
                                                   loopsamples=20,
                                                   resultfile='target')
            soup.yasara.Experiment("On")
            soup.yasara.Wait("Expend")

            if os.path.isfile(error_path):
                with open(error_path, 'r') as f:
                    raise ModelRunError(f.read())
            elif os.path.isfile(error_scene_path):
                raise ModelRunError("yasara exited with an error")

            if not os.path.isfile(output_yob_path):
                raise ModelRunError("yasara generated no output yob")

            model_path = os.path.join(work_dir_path, 'target.pdb')
            soup.yasara.SavePDB(soup.template_obj, model_path)

            self._write_selected_targets(chain_alignments, os.path.join(work_dir_path, 'selected-targets.txt'))

            tar_path = model_storage.get_tar_path(soup.get_main_target_sequence(),
                                                  soup.target_species_id,
                                                  main_domain_alignment,
                                                  TemplateID(soup.template_pdbid,
                                                             soup.main_target_chain_id))
            with tarfile.open(tar_path, mode="w:gz") as ar:
                ar.add(work_dir_path, arcname=model_name)

            return tar_path
        except RuntimeError as e:
            self._log_additional_error_info(e, chain_alignments)

            if os.path.isfile(error_path):
                with open(error_path, 'r') as f:
                    raise ModelRunError(f.read())
            elif os.path.isfile(error_scene_path):
                raise ModelRunError("yasara exited with an error")
            else:
                raise e
        finally:
            if os.path.isdir(work_dir_path):
                shutil.rmtree(work_dir_path)

    def _write_selected_targets(self, alignments_per_chain, path):
        with open(path, 'w') as f:
            for chain_id in alignments_per_chain:
                if alignments_per_chain[chain_id].target_id is not None:
                    f.write("%s: %s\n" % (chain_id, alignments_per_chain[chain_id].target_id))

    def _log_additional_error_info(self, err, alignments):
        message = str(err)
        if P_TEMPLATE_DIFFERENT.search(message):

            for chain_id in soup.get_chain_ids():
                sequence = soup.get_sequence(chain_id)
                _log.error("template chain {} has sequence {}".format(chain_id, sequence))
                _log.error("vs sequencemol {}".format(soup.yasara.SequenceMol("obj %i and mol %s and protein" % (soup.template_obj, chain_id))))


P_TEMPLATE_DIFFERENT = re.compile(r"The template object [0-9]+ does not contain all the residues specified in the sequence alignment")

modeler = Modeler()
