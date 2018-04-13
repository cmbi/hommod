from hommod.models.aminoacid import AminoAcid


def parse_seqres_from_string(s):
    sequences = {}
    for line in s.split('\n'):
        if line.startswith('SEQRES'):
            chain_id = line[11]
            if chain_id not in sequences:
                sequences[chain_id] = []

            for code in line[19:].split():
                sequences[chain_id].append(AminoAcid.from_three_letter_code(code))

    return sequences
