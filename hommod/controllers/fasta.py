


def write_fasta(path, d):
    with open(path, 'w') as f:
        for key in d:
            f.write('>%s\n%s\n' % (key, d[key]))


def parse_fasta_from_string(s):
    d = {}
    current_id = None
    for line in s.split('\n'):
        if line.startswith('>'):
            current_id = line[1:].strip()
            d[current_id] = ''
        else:
            d[current_id] += line.strip()
    return d


def parse_fasta(path):
    d = {}
    current_id = None
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                current_id = line[1:].strip()
                d[current_id] = ''
            else:
                d[current_id] += line.strip()
    return d
