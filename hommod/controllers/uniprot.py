from hommod.controllers.fasta import FastaIterator
from hommod.models.error import InitError


class Uniprot:
    def __init__(self, fasta_paths=None):
        self.fasta_paths = fasta_paths

    def get_sequence(self, ac):
        if self.fasta_paths is None:
            raise InitError("fasta paths not set")

        for fasta_path in self.fasta_paths:
            with FastaIterator(fasta_path) as fasta:
                for id_, seq in fasta:
                    if id_.split('|')[1] == ac:
                        return seq

        raise ValueError("sequence not found in uniprot: {}".format(ac))

uniprot = Uniprot()
