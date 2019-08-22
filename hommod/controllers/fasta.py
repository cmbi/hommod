


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


class FastaIterator:
    def __init__(self, filepath):
        self.filepath = filepath

    def __enter__(self):
        self._file = open(self.filepath, 'r')
        return self

    def __exit__(self, type_, value, traceback):
        self._file.close()

    def __iter__(self):
        return self

    def __next__(self):
        id_ = None
        sequence = ""
        while True:
            line = self._file.readline()

            if len(line) <= 0:
                if id_ is not None:
                    return id_, sequence
                else:
                    raise StopIteration

            elif line.startswith('>'):
                if id_ is not None:
                    self._file.seek(self._file.tell() - len(line), 0)
                    return id_, sequence

                id_ = line[1:].split()[0]

            elif id_ is not None:
                sequence += line

        return id_, sequence
