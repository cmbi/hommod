


class TemplateID:
    def __init__(self, pdbid, chain_id):
        self.pdbid = pdbid.lower()
        self.chain_id = chain_id

    def __eq__(self, other):
        return self.pdbid.lower() == other.pdbid.lower() and \
               self.chain_id == other.chain_id

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.pdbid.lower(), self.chain_id))

    def __repr__(self):
        return "%s-%s" % (self.pdbid, self.chain_id)
