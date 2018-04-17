import os

from hommod.models.error import InitError


class BlackLister:
    def __init__(self, file_path=None):
        self.file_path = file_path

    def is_blacklisted(self, pdbid):
        if self.file_path is None:
            raise InitError("blacklist file not set")

        if os.path.isfile(self.file_path):
            with open(self.file_path, 'r') as f:
                list_ = f.read().split()
                return pdbid in list_

        return False

    def blacklist(self, pdbid):
        if self.file_path is None:
            raise InitError("blacklist file not set")

        if not is_blacklisted(pdbid):
            with open(self.file_path, 'a') as f:
                f.write('%s\n' % pdbid)

blacklister = BlackLister()
