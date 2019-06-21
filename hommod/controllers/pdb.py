from hommod.models.aminoacid import AminoAcid


def parse_seqres_from_string(s):
    sequences = {}
    for line in s.split('\n'):
        if line.startswith('SEQRES'):
            chain_id = line[11]
            if chain_id not in sequences:
                sequences[chain_id] = []

            for code in line[19:].split():
                if len(code) == 3 and code != 'HOH':
                    sequences[chain_id].append(AminoAcid.from_three_letter_code(code))

    return sequences


def parse_chain_order_from_string(s):
    order = []
    for line in s.split('\n'):
        if line.startswith('ATOM'):
            code = line[17: 20].strip()
            if AminoAcid.from_three_letter_code(code).name != "Other":
                chain_id = line[21]

                if len(order) == 0 or chain_id not in order:
                    order.append(chain_id)

    return order
