import re

from hommod.models.range import SequenceRange
from hommod.controllers.sequence import is_amino_acid_char


class Alignment:
    def __init__(self, aligned_sequences):

        # Alignments are assumed to have '-' for gaps.
        self.aligned_sequences = aligned_sequences

    def __repr__(self):
        s = "\n"
        n = 100
        m = 20
        l = len(self.aligned_sequences.values()[0])
        for i in range(0, l, n):
            f = min(i + n, l)
            for k in self.aligned_sequences:
                s += (m - len(k)) * ' ' + k + ': ' + self.aligned_sequences[k][i: f] + '\n'
            s += '\n'

        return s

    def get_sequence(self, id_):
        return self.aligned_sequences[id_].replace('-','')

    def count_aligned_residues(self, id1, id2):
        nalign = 0
        for i in range(len(self.aligned_sequences[id1])):
            if is_amino_acid_char(self.aligned_sequences[id1][i]) and \
                    is_amino_acid_char(self.aligned_sequences[id2][i]):
                nalign += 1
        return nalign

    def get_percentage_identity(self, id1, id2):
        nalign = 0
        nid = 0
        for i in range(len(self.aligned_sequences[id1])):
            if is_amino_acid_char(self.aligned_sequences[id1][i]) and \
                    is_amino_acid_char(self.aligned_sequences[id2][i]):
                nalign += 1
                if self.aligned_sequences[id1][i] == self.aligned_sequences[id2][i]:
                    nid += 1
        if nalign > 0:
            return (100.0 * nid) / nalign
        else:
            return 0.0

    def as_dict(self):
        return {key: self.aligned_sequences[key] for key in self.aligned_sequences}


pdbid_p = re.compile(r'[0-9][0-9a-zA-Z]{3}_[0-9a-zA-Z]')

class BlastAlignment(Alignment):
    def __init__(self, hit_id, full_query_sequence,  databank,
                 query_start, query_end, query_alignment,
                 subject_start, subject_end, subject_alignment):
        Alignment.__init__(self, {'query': query_alignment,
                                  'subject': subject_alignment})

        self.hit_id = hit_id
        self.full_query_sequence = full_query_sequence
        self.databank = databank

        # Start and end are inclusive, residue numbers start from 1.
        self.query_start = query_start
        self.query_end = query_end
        self.subject_start = subject_start
        self.subject_end = subject_end

    def __repr__(self):
        return self.hit_id + '\n' + Alignment.__repr__(self)

    def get_query_range(self):
        return SequenceRange(self.query_start - 1, self.query_end, self.full_query_sequence)

    def get_hit_type(self):
        if pdbid_p.match(self.hit_id):
            return 'pdb'
        elif '|' in self.hit_id:
            s = self.hit_id.split('|')
            return s[0]
        else:
            raise ValueError("Cannot parse hit id {}".format(self.hit_id))

    def get_hit_accession_code(self):
        if pdbid_p.match(self.hit_id):
            return self.hit_id.split('_')[0]
        elif '|' in self.hit_id:
            s = self.hit_id.split('|')
            return s[1]
        else:
            raise ValueError("Cannot parse hit id {}".format(self.hit_id))

    def get_hit_chain_id(self):
        if pdbid_p.match(self.hit_id):
            return self.hit_id.split('_')[1]
        elif '|' in self.hit_id:
            s = self.hit_id.split('|')
            if s[0] not in ['pdb', 'pdbfinder']:
                raise ValueError("Not a pdb hit id: {}".format(self.hit_id))

            return s[2]
        else:
            raise ValueError("Cannot parse hit id {}".format(self.hit_id))

    @property
    def query_alignment(self):
        return self.aligned_sequences['query']

    @query_alignment.setter
    def query_alignment(self, value):
        self.aligned_sequences['query'] = value

    @property
    def subject_alignment(self):
        return self.aligned_sequences['subject']

    @query_alignment.setter
    def subject_alignment(self, value):
        self.aligned_sequences['subject'] = value

    def get_percentage_identity(self):
        return Alignment.get_percentage_identity(self, 'query', 'subject')

    def count_aligned_residues(self):
        return Alignment.count_aligned_residues(self, 'query', 'subject')

    def get_percentage_coverage(self):
        nalign = 0
        nq = 0
        for i in range(len(self.aligned_sequences['query'])):
            if is_amino_acid_char(self.aligned_sequences['query'][i]):
                nq += 1
                if is_amino_acid_char(self.aligned_sequences['subject'][i]):
                    nalign += 1
        return (100.0 * nalign) / nq

    def is_query_residue_covered(self, residue_number):
        n = self.query_start

        for i in range(len(self.aligned_sequences['query'])):
            if is_amino_acid_char(self.aligned_sequences['query'][i]):
                if n == residue_number and \
                        is_amino_acid_char(self.aligned_sequences['subject'][i]):
                    return True
                n += 1

        return False


class TargetTemplateAlignment(Alignment):
    def __init__(self, target_alignment, template_alignment):
        Alignment.__init__(self, {'target': target_alignment,
                                  'template': template_alignment})
        # These may be set later:
        self.target_id = None
        self.template_id = None

    @property
    def target_alignment(self):
        return self.aligned_sequences['target']

    @target_alignment.setter
    def target_alignment(self, value):
        self.aligned_sequences['target'] = value

    @property
    def template_alignment(self):
        return self.aligned_sequences['template']

    @template_alignment.setter
    def template_alignment(self, value):
        self.aligned_sequences['template'] = value

    def get_target_sequence(self):
        return self.get_sequence('target')

    def get_template_sequence(self):
        return self.get_sequence('template')

    def get_percentage_identity(self):
        return Alignment.get_percentage_identity(self, 'target', 'template')

    def count_aligned_residues(self):
        return Alignment.count_aligned_residues(self, 'target', 'template')

    def get_percentage_coverage(self):
        nalign = 0
        ntar = 0
        for i in range(len(self.aligned_sequences['target'])):
            if is_amino_acid_char(self.aligned_sequences['target'][i]):
                ntar += 1
                if is_amino_acid_char(self.aligned_sequences['template'][i]):
                    nalign += 1
        return (100.0 * nalign) / ntar

    def get_covered_template_residues_indices(self):
        n = 0
        covered = []
        for i in range(len(self.aligned_sequences['target'])):
            if is_amino_acid_char(self.aligned_sequences['template'][i]):
                if is_amino_acid_char(self.aligned_sequences['target'][i]):
                    covered.append(n)
                n += 1
        return covered

    def is_target_residue_covered(self, residue_number):
        n = 1
        for i in range(len(self.aligned_sequences['target'])):
            if is_amino_acid_char(self.aligned_sequences['target'][i]):
                if n == residue_number and \
                        is_amino_acid_char(self.aligned_sequences['template'][i]):
                    return True
                n += 1

        return False

    def get_relative_span(self):
        """
        Tells the starting position of 'target' relative to
        the starting position of 'template'.
        """

        i = 0 
        while not is_amino_acid_char(self.target_alignment[i]):
            i += 1

        start = len(self.template_alignment[:i].replace('-', ''))

        i = len(self.target_alignment)
        while not is_amino_acid_char(self.target_alignment[i - 1]):
            i -= 1

        end = len(self.template_alignment[:i].replace('-', ''))

        return SequenceRange(start, end, self.get_template_sequence())

class DomainAlignment(TargetTemplateAlignment):
    def __init__(self, target_alignment, template_alignment, range_, template_id):
        TargetTemplateAlignment.__init__(self, target_alignment, template_alignment)
        self.range = range_
        self.template_id = template_id

    def __repr__(self):
        return "{} {}\n".format(self.template_id, self.range) + TargetTemplateAlignment.__repr__(self)
