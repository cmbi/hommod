

def is_protein_sequence(s):
    return all([is_amino_acid_char(s[i]) for i in range(len(s))])

def is_nucleic_acid_sequence(s):
    return all([is_nucleotide_char(s[i]) for i in range(len(s))])

def is_amino_acid_char(c):
    return c in ['A','B','C','D','E','F','G','H','I','J','K',
                 'L','M','N','O','P','Q','R','S','T','U','V',
                 'W','X','Y','Z']

def is_nucleotide_char(c):
    return c in ['a', 't', 'c', 'u', 'g']
