import sys
import imp
import os

from hommod.models.error import ModelRunError, InitError
from hommod.models.aminoacid import AminoAcid
from hommod.models.residue import ModelingResidue
from hommod.controllers.yasara import YasaraContext


class ModelingContext:

    def __init__(self, yasara_dir):
        self.yasara = YasaraContext(yasara_dir)

        self.template_obj = None
        self.template_pdbid = None
        self.main_target_chain_id = None
        self.target_species_id = None
        self.target_sequences = {}

    def set_main_target(self, main_target_sequence, target_species_id, main_target_chain_id):
        self.target_species_id = target_species_id
        self.main_target_chain_id = main_target_chain_id
        self.target_sequences[self.main_target_chain_id] = main_target_sequence

    def get_main_target_sequence(self):
        return self.target_sequences[self.main_target_chain_id]

    def get_chain_ids(self):
        if self.template_obj is None:
            raise ModelRunError("template object is not set")

        chain_ids = self.yasara.ListMol('obj %i and protein' % self.template_obj, 'MOL')
        return chain_ids

    def delete_chain(self, chain_id):
        if self.template_obj is None:
            raise ModelRunError("template object is not set")

        self.yasara.DelMol('obj %i and protein and mol %s' % (self.template_obj, chain_id))

    def get_sequence(self, chain_id):
        sequence = ""
        for residue in self.get_residues(chain_id):
            sequence += residue.amino_acid.letter

        return sequence

    def get_residues(self, chain_id):
        if self.template_obj is None:
            raise ModelRunError("template object is not set")

        residues = []
        for s in self.yasara.ListAtom("obj %i and mol %s and protein" % (self.template_obj, chain_id),
                                      "RESNUM RESNAME ATOMNAME ATOMNUM"):
            resnum, resname, atomname, atomnum = s.split()

            if len(residues) <= 0 or residues[-1].residue_number != resnum:
                amino_acid = AminoAcid.from_three_letter_code(resname)
                residues.append(ModelingResidue(resnum, amino_acid))

            residues[-1].atom_numbers[atomname] = int(atomnum)

        return residues

    def list_interacting_chains(self, chain_id):
        return self.yasara.ListMol(
                 'protein and obj %i and not mol %s with distance<4.5 from obj %i mol %s' %
                 (self.template_obj, chain_id, self.template_obj, chain_id),
                 'MOL')

    def residue_interacts_with(self, residue, interacting_residues):
        if self.template_obj is None:
            raise ModelRunError("template object is not set")

        if 'CA' not in residue.atom_numbers:
            return False

        interacting_residues = list(filter(lambda res: 'CA' in res.atom_numbers, interacting_residues))

        atoms = self.yasara.ListAtom("CA and atom %i-%i and obj %i with distance<6 from %i" %
                                     (interacting_residues[0].atom_numbers['CA'],
                                      interacting_residues[-1].atom_numbers['CA'],
                                      self.template_obj, residue.atom_numbers['CA']))
        return len(atoms) > 0

    def get_secondary_structure(self, chain_id):
        if self.template_obj is None:
            raise ModelRunError("template object is not set")

        return ''.join(self.yasara.SecStrRes('obj %i and protein and mol %s' % (self.template_obj, chain_id)))


