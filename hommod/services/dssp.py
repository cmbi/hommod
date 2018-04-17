import os
import logging


_log = logging.getLogger(__name__)


class DsspService:
    def __init__(self, dssp_dir=None):
        self.dssp_dir = dssp_dir

    def has_secondary_structure(self, template_id):
        try:
            dssp_str = self._get_dssp(template_id.pdbid)
        except:
            return False

        data = self._parse_dssp(dssp_str)

        return template_id.chain_id in data

    def get_sequence(self, template_id):
        dssp_str = self._get_dssp(template_id.pdbid)
        data = self._parse_dssp(dssp_str)

        return data[template_id.chain_id][0]

    def get_secondary_structure(self, template_id):
        dssp_str = self._get_dssp(template_id.pdbid)
        data = self._parse_dssp(dssp_str)

        return data[template_id.chain_id][1]

    def _get_dssp(self, pdbid):
        file_path = os.path.join(self.dssp_dir, '%s.dssp' % pdbid.lower())
        with open(file_path, 'r') as f:
            return f.read()

    def _parse_dssp(self, dssp_str):
        data = {}
        for line in dssp_str.split('\n'):

            stripped_line = line.strip()
            if stripped_line.endswith('.') or stripped_line.startswith('#') or len(stripped_line) <= 0:
                continue

            chain_id = line[11]
            amino_acid = line[13]
            secstr = line[16]

            if not chain_id.isalpha():
                continue

            if chain_id not in data:
                data[chain_id] = ['', '', {}] 

            i = len(data[chain_id][0])

            if amino_acid.islower():  # disulfid bridge

                if amino_acid not in data[chain_id][2]:
                    data[chain_id][2][amino_acid] = []
                data[chain_id][2][amino_acid].append(i)

                amino_acid = 'C' 

            data[chain_id][0] += amino_acid
            data[chain_id][1] += secstr

        return data

dssp = DsspService()
