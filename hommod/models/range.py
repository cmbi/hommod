


class SequenceRange:
    def __init__(self, start, end, sequence):
        if start >= end:
            raise ValueError("start must be smaller than end")

        self.sequence = sequence
        self.start = start
        self.end = end
        self.ac = None

    def is_left_from(self, other):

        if self.start < other.start:
            return True
        elif self.start == other.start:
            return self.end < other.end
        else:
            return False

    def is_right_from(self, other):

        if self.end > other.end:
            return True
        elif self.end == other.end:
            return self.start > other.start
        else:
            return False

    def get_distance_from(self, other):
        if self.start < other.start:
            if self.end < other.start:
                return other.start - self.end
            else:  # overlap
                return 0

        else:  # other.start < self.start

            if other.end < self.start:
                return self.start - other.end
            else:  # overlap
                return 0

    def get_intersection(self, other):
        if self.sequence != other.sequence:
            raise ValueError("Not from the same Sequence")

        return SequenceRange(max(self.start, other.start),
                             min(self.end, other.end),
                             self.sequence)


    def __sub__(self, value):
        return SequenceRange(self.start - value, self.end - value, self.sequence)

    def __eq__(self, other):
        return (self.start == other.start and self.end == other.end and self.sequence == other.sequence)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.start, self.end, self.sequence))

    def __repr__(self):
        return '%i-%i' % (self.start, self.end)

    def get_length(self):
        return self.end - self.start

    def get_sub_sequence(self):
        return self.sequence[self.start: self.end]

    def includes_residue(self, resnum):
        residue_index = resnum - 1
        return residue_index >= self.start and residue_index < self.end

    def get_percentage_overlap(self, other):
        if not self.overlaps_with(other):
            return 0.0

        count_overlap = min(self.end, other.end) - max(self.start, other.start)
        if count_overlap <= 0:
            return 0.0

        return  (100.0 * count_overlap) / min(self.get_length(), other.get_length())

    def overlaps_with(self, other):
        if self.sequence != other.sequence:
            return False

        # given: start < end
        return (self.end >= other.start and self.start <= other.end or
                other.end >= self.start and other.start <= self.end)

    def encloses(self, other):
        if self.sequence != other.sequence:
            return False

        # all of 'other' lies in 'self'
        return (self.start <= other.start and self.end >= other.end)

    def merge_with(self, other):
        if not self.overlaps_with(other):
            raise ValueError("not overlapping")

        return SequenceRange(min(self.start, other.start),
                             max(self.end, other.end),
                             self.sequence)
