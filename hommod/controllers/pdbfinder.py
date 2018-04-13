class ResultNode:
    def __init__(self, key, value):
        self.key = key
        self.value = value
        self.children = []


    def find(self, key):
        return list(filter(lambda child: child.key == key, self.children))

    def find_one(self, key):
        for child in self.children:
            if child.key == key:
                return child
        return None

    def __repr__(self):
        return "{}: {}; with children {}".format(self.key, self.value, self.children)


class _NodeBuilder:
    def __init__(self, key, value):
        self.key = key
        self.value = value

        self._nodes = []
        self._last_indent = 0

    def add(self, indent, key, value):

        if indent > self._last_indent and len(self._nodes) > 0:
            self._nodes[-1].add(indent, key, value)
        else:
            self._nodes.append(_NodeBuilder(key, value))
            self._last_indent = indent

    def get_result(self):
        r = ResultNode(self.key, self.value)
        r.children = [n.get_result() for n in self._nodes]
        return r


def _get_indent(line):
    i = 0
    while line[i].isspace():
        i += 1

    return i


def parse_pdbfinder(file_path):

    entries = []

    current_entry = None

    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('//'):  # comment
                continue

            i = line.find(':')
            if i < 0:  # no colon in line
                continue

            indent = _get_indent(line)
            key = line[:i].strip()
            value = line[i + 1:].strip()

            if key == 'ID' and indent == 0:
                if current_entry is not None:
                    entries.append(current_entry.get_result())
                current_entry = _NodeBuilder(None, None)

            current_entry.add(indent, key, value)

        if current_entry is not None:
            entries.append(current_entry.get_result())

    return entries
