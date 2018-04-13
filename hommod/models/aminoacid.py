import logging


_log = logging.getLogger(__name__)

class AminoAcid:
    @staticmethod
    def from_three_letter_code(code):
        matches = [aa for aa in _AMINO_ACIDS if aa.code.upper() == code.upper()]
        if len(matches) <= 0:
            return AminoAcid.from_three_letter_code('OTH')
        assert len(matches) == 1
        return matches[0]

    @staticmethod
    def from_one_letter_code(code):
        matches = [aa for aa in _AMINO_ACIDS if aa.code.upper() == code.upper()]
        if len(matches) <= 0:
            return AminoAcid.from_one_letter_code('X')
        assert len(matches) == 1
        return matches[0]

    def __init__(self, code, letter, name):
        self._code = code
        self._letter = letter
        self._name = name

    @property
    def code(self):
        return self._code

    @property
    def letter(self):
        return self._letter

    @property
    def name(self):
        return self._name


_AMINO_ACIDS = (
    AminoAcid('ALA', 'A', 'Alanine'),
    AminoAcid('ARG', 'R', 'Arginine'),
    AminoAcid('ASN', 'N', 'Asparagine'),
    AminoAcid('ASP', 'D', 'Aspartic Acid'),
    AminoAcid('CYS', 'C', 'Cysteine'),
    AminoAcid('GLU', 'E', 'Glutamic Acid'),
    AminoAcid('GLN', 'Q', 'Glutamine'),
    AminoAcid('GLY', 'G', 'Glycine'),
    AminoAcid('HIS', 'H', 'Histidine'),
    AminoAcid('ILE', 'I', 'Isoleucine'),
    AminoAcid('LEU', 'L', 'Leucine'),
    AminoAcid('LYS', 'K', 'Lysine'),
    AminoAcid('MET', 'M', 'Methionine'),
    AminoAcid('PHE', 'F', 'Phenylalanine'),
    AminoAcid('PRO', 'P', 'Proline'),
    AminoAcid('SER', 'S', 'Serine'),
    AminoAcid('THR', 'T', 'Threonine'),
    AminoAcid('TRP', 'W', 'Tryptophan'),
    AminoAcid('TYR', 'Y', 'Tyrosine'),
    AminoAcid('VAL', 'V', 'Valine'),

    # Non standard
    AminoAcid('SEC', 'U', 'Selenocysteine'),
    AminoAcid('PYL', 'O', 'Pyrrolysine'),

    # Unknown
    AminoAcid('OTH', 'X', 'Other')
)
