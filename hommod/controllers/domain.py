import logging
import traceback

from hommod.controllers.rost import get_min_identity
from hommod.controllers.blast import blaster
from hommod.controllers.blacklist import blacklister
from hommod.services.interpro import interpro
from hommod.services.dssp import dssp
from hommod.models.template import TemplateID
from hommod.models.range import SequenceRange
from hommod.models.align import DomainAlignment
from hommod.models.error import InitError
from hommod.controllers.sequence import is_amino_acid_char
from hommod.controllers.kmad import kmad_aligner
from hommod.controllers.log import ModelLogger


_log = logging.getLogger(__name__)


class DomainAligner:
    def __init__(self, forbidden_interpro_domains=None,
                 similar_ranges_min_overlap_percentage=None,
                 similar_ranges_max_length_difference_percentage=None,
                 template_blast_databank=None,
                 min_percentage_coverage=None,
                 highly_homologous_percentage_identity=None):

            self.forbidden_interpro_domains = forbidden_interpro_domains
            self.similar_ranges_min_overlap_percentage = similar_ranges_min_overlap_percentage
            self.similar_ranges_max_length_difference_percentage = similar_ranges_max_length_difference_percentage
            self.template_blast_databank = template_blast_databank
            self.min_percentage_coverage = min_percentage_coverage
            self.highly_homologous_percentage_identity = highly_homologous_percentage_identity

    def get_domain_alignments(self, target_sequence, require_resnum=None, template_id=None):

        ModelLogger.get_current().add("getting domain alignments for sequence {}, resnum {}, template {}"
                                      .format(target_sequence, require_resnum, template_id))

        if self.min_percentage_coverage is None:
            raise InitError("min percentage coverage is not set")

        interpro_ranges = interpro.get_domain_ranges(target_sequence)
        _log.debug("{} ranges from interpro".format(len(interpro_ranges)))

        sample_ranges = self._filter_forbidden_ranges(interpro_ranges)

        if require_resnum is not None:
            sample_ranges = list(filter(lambda r: r.includes_residue(require_resnum), sample_ranges))
            _log.debug("{} ranges have residue {}".format(len(sample_ranges), require_resnum))

        # Add the whole sequence as a range too:
        sample_ranges.append(SequenceRange(0, len(target_sequence), target_sequence))

        ok_ranges_alignments = {}
        best_ranges_alignments = {}
        checked_ranges = []

        while len(sample_ranges) > 0:

            merged_sample_ranges = self._merge_similar_ranges(sample_ranges)

            _log.debug("sampling {} ranges".format(len(merged_sample_ranges)))

            # Check the largest ranges first. If that yields, then the smaller ones don't matter.
            for range_ in sorted(merged_sample_ranges, key=lambda r: r.get_length(), reverse=True):

                if range_ in checked_ranges:
                    continue  # already passed this one
                checked_ranges.append(range_)

                if any([r.encloses(range_) for r in best_ranges_alignments]):
                    continue  # we already have a larger enclosing range

                # These can differ per range:
                best_hit = None
                last_resort_hit = None

                ModelLogger.get_current().add("examining range {}".format(range_))

                hit_candidates = self._get_hits(range_, template_id)

                _log.debug('trying range: {} against {} hits'.format(range_, len(hit_candidates)))

                for hit_candidate in hit_candidates:
                    hit_range = hit_candidate.get_query_range()
                    if require_resnum is not None:
                        if not hit_candidate.is_query_residue_covered(require_resnum):
                            _log.debug("hit with {} on {} does not cover residue {}"
                                       .format(hit_candidate.get_hit_accession_code(),
                                               hit_range, require_resnum))
                            continue

                    if self._alignment_ok_for_range(range_, hit_candidate):
                        _log.debug("hit with {} {} is ok"
                                   .format(hit_candidate.get_hit_accession_code(), hit_range))

                        # This range made an OK alignment, so at least store it for later usage:
                        template_id = TemplateID(hit_candidate.get_hit_accession_code(),
                                                 hit_candidate.get_hit_chain_id())
                        ok_ranges_alignments[hit_range] = DomainAlignment(hit_candidate.query_alignment,
                                                                          hit_candidate.subject_alignment,
                                                                          hit_range, template_id)

                        ModelLogger.get_current().add("found a hit with {} covering range {}:\n{}"
                                                      .format(template_id, hit_range, hit_candidate))


                        if hit_candidate.get_percentage_coverage() > self.min_percentage_coverage:

                            _log.debug("coverage is high enough for {} {}"
                                       .format(hit_candidate.get_hit_accession_code(), hit_range))

                            if best_hit is None or self._is_better_than(hit_candidate, best_hit):

                                _log.debug("{} is better than {}".format(hit_candidate, best_hit))
                                ModelLogger.get_current().add("{} is better than {}".format(hit_candidate, best_hit))

                                best_hit = hit_candidate
                        else:
                            last_resort_hit = hit_candidate

                if best_hit is None:
                    best_hit = last_resort_hit

                if best_hit is not None:

                    # Remove any smaller ranges that this one encloses:
                    best_ranges_alignments = self._remove_enclosing(range_, best_ranges_alignments)

                    template_id = TemplateID(best_hit.get_hit_accession_code(),
                                             best_hit.get_hit_chain_id())

                    hit_range = best_hit.get_query_range()
                    _log.debug("passing best hit with template {} with range {}".format(template_id, hit_range))

                    best_ranges_alignments[hit_range] = DomainAlignment(best_hit.query_alignment,
                                                                        best_hit.subject_alignment,
                                                                        hit_range, template_id)
                else:
                    _log.debug("no hit for range {}".format(range_))

            # After iterating the sample ranges, prepare for the next round:
            sample_ranges = self._clean_search_space(checked_ranges, sample_ranges, ok_ranges_alignments)

        return best_ranges_alignments.values()

    def _remove_enclosing(self, range_, ranges_alignments):
        new_ranges_alignments = {}
        for r in ranges_alignments:
            if not range_.encloses(r):
                new_ranges_alignments[r] = ranges_alignments[r]
        return new_ranges_alignments

    def _alignment_ok_for_range(self, range_, alignment):

        pid = alignment.get_percentage_identity()
        nalign = alignment.count_aligned_residues()
        pcover = (100.0 * nalign) / range_.get_length()

        highly_homolgous = pid >= self.highly_homologous_percentage_identity and \
                           range_.get_length() == len(range_.sequence)

        _log.debug("alignment with {}: pid={}, nalign={}, pcover={}"
                   .format(range_, pid, nalign, pcover))

        return pid >= get_min_identity(nalign) and \
               (pcover >= self.min_percentage_coverage or highly_homolgous)

    def _find_shared_hits_ranges(self, ranges_alignments):
        """
        If two ranges have the same template, they are put together
        in the returned dictionary as a list.
        """
        # sort by template id
        d = {}
        for r in ranges_alignments:
            if ranges_alignments[r].template_id not in d:
                d[ranges_alignments[r].template_id] = []
            d[ranges_alignments[r].template_id].append(r)
 
        # Forget those blast hits that only occur once:
        dr = {}
        for template_id in d.keys():
            if len(d[template_id]) > 1:
                dr[template_id] = d[template_id]

        return dr

    def _clean_search_space(self, checked_ranges, sample_ranges, ok_ranges_alignments):

        # See if we can merge ranges that have
        # the same template in their blast hits:
        checked_ranges = self._remove_duplicate_ranges(checked_ranges + sample_ranges)
        sample_ranges = []
        shared_hits_ranges = self._find_shared_hits_ranges(ok_ranges_alignments)
        for template_id in shared_hits_ranges:

            ranges = shared_hits_ranges[template_id]

            for i in range(len(ranges)):
                overlapping_indices = []
                for j in range(len(ranges)):
                    if j != i and ranges[j].overlaps_with(ranges[i]):
                        overlapping_indices.append(j)

                for j in overlapping_indices:
                    percentage_overlap = ranges[i].get_percentage_overlap(ranges[j])
                    percentage_length_difference = 100.0 * (abs(ranges[i].get_length() - ranges[j].get_length()) /
                                                            max(ranges[i].get_length(), ranges[j].get_length()))

                    merged = ranges[i].merge_with(ranges[j])

                    # Merge only if:
                    # - the ranges are close together
                    # - the merge has not already been done
                    # - the intersecting parts of the ranges align to
                    #   the template in exactly the same way
                    if merged not in checked_ranges:

                        alignment_i = ok_ranges_alignments[ranges[i]]
                        alignment_j = ok_ranges_alignments[ranges[j]]
                        template_secstr = dssp.get_secondary_structure(template_id)
                        template_sequence = dssp.get_sequence(template_id)
                        try:
                            alignment_m = kmad_aligner.align(template_sequence, template_secstr,
                                                             merged.get_sub_sequence())
                        except:
                            _log.warn(traceback.format_exc())

                            # If kmad fails, then skip this one :(
                            continue

                        intersected = ranges[i].get_intersection(ranges[j])

                        intersect_template_sequence_i = \
                            self._get_template_sequence_in_target_range(alignment_i, intersected - ranges[i].start)

                        intersect_template_sequence_j = \
                            self._get_template_sequence_in_target_range(alignment_j, intersected - ranges[j].start)

                        intersect_template_sequence_m = \
                            self._get_template_sequence_in_target_range(alignment_m, intersected - merged.start)

                        if intersect_template_sequence_i == intersect_template_sequence_m and \
                                intersect_template_sequence_j == intersect_template_sequence_m:
                            sample_ranges.append(merged)

        return sample_ranges

    def _get_template_sequence_in_target_range(self, alignment, target_range):

        """
        For the given target sequence range, this function returns the
        corresponding sequence of the template, according to the given alignment.
        """

        target_start = 0
        while not is_amino_acid_char(alignment.target_alignment[target_start]):
            target_start += 1

        start = target_start
        n_aa = 0
        while n_aa < target_range.start and start < len(alignment.target_alignment):
            if is_amino_acid_char(alignment.target_alignment[start]):
                n_aa += 1
            start += 1

        end = start
        n_aa = 0
        while n_aa < target_range.get_length() and end < len(alignment.target_alignment):
            if is_amino_acid_char(alignment.target_alignment[end]):
                n_aa += 1
            end += 1

        return alignment.template_alignment[start: end]

    def _get_hits(self, range_, template_id):
        if self.template_blast_databank is None:
            raise InitError("blast databank is not set")

        blast_hits = blaster.blastp(range_.get_sub_sequence(), self.template_blast_databank)
        _log.debug("{} blast hits to filter".format(len(blast_hits)))

        good_hits = []
        for hit_id in blast_hits:
            for alignment in blast_hits[hit_id]:
                hit_template_id = TemplateID(alignment.get_hit_accession_code(),
                                             alignment.get_hit_chain_id())
                if template_id is not None and hit_template_id != template_id:
                    continue

                if template_id is None and blacklister.is_blacklisted(alignment.get_hit_accession_code()):
                    continue

                if not dssp.has_secondary_structure(hit_template_id):
                    continue

                # Replace the blast hit's alignment with the kmad alignment.
                template_secstr = dssp.get_secondary_structure(hit_template_id)
                template_sequence = dssp.get_sequence(hit_template_id)
                try:
                    kmad_alignment = kmad_aligner.align(template_sequence, template_secstr,
                                                        range_.get_sub_sequence())
                except:
                    _log.warn(traceback.format_exc())

                    # If kmad fails, then skip this one :(
                    continue
                alignment.full_query_sequence = range_.sequence
                alignment.query_start = range_.start + 1
                alignment.query_end = range_.end
                alignment.subject_start = 1
                alignment.subject_end = len(template_sequence)
                alignment.query_alignment = kmad_alignment.target_alignment
                alignment.subject_alignment = kmad_alignment.template_alignment

                if alignment.get_percentage_identity() >= get_min_identity(alignment.count_aligned_residues()):
                    good_hits.append(alignment)

        return good_hits

    def _is_better_than(self, hit, other_hit):
        _log.debug("compare new {} {}% with current best {} {}%"
                   .format(hit.hit_id, hit.get_percentage_identity(),
                           other_hit.hit_id, other_hit.get_percentage_identity()))

        return hit.get_percentage_identity() >= other_hit.get_percentage_identity() and \
               hit.count_aligned_residues() >= other_hit.count_aligned_residues()

    def _merge_similar_ranges(self, ranges):
        if self.similar_ranges_min_overlap_percentage is None or \
                self.similar_ranges_max_length_difference_percentage is None:
            raise InitError("similar range percentages not set")

        ranges = sorted(ranges, key=lambda r: r.start)

        i = 0
        while i < len(ranges):
            overlapping_indices = []
            for j in range(i + 1, len(ranges)):
                if ranges[j].overlaps_with(ranges[i]):
                    overlapping_indices.append(j)

            # important, rightmost must go first!
            # Because we're going to remove ranges from the list.
            overlapping_indices = sorted(overlapping_indices, reverse=True)
            for j in overlapping_indices:

                percentage_overlap = ranges[i].get_percentage_overlap(ranges[j])
                percentage_length_difference = ((100.0 * abs(ranges[i].get_length() - ranges[j].get_length())) /
                                                max(ranges[i].get_length(), ranges[j].get_length()))

                if percentage_overlap > self.similar_ranges_min_overlap_percentage and \
                        percentage_length_difference < self.similar_ranges_max_length_difference_percentage:

                    # Replace the two ranges by a single merged one:
                    _log.debug("merging {} with {}, they have {} % length difference"
                               .format(ranges[i], ranges[j], percentage_length_difference))
                    merged = ranges[i].merge_with(ranges[j])
                    ranges = (ranges[:i] + [merged] + ranges[i + 1: j] + ranges[j + 1:])
            i += 1

            # Make list shorter to save time:
            ranges = self._remove_duplicate_ranges(ranges)

        return ranges

    def _remove_duplicate_ranges(self, ranges):
        i = 0
        while i < len(ranges):
            j = i + 1
            while j < len(ranges):
                if ranges[i] == ranges[j]:
                    if (j + 1) < len(ranges):
                        ranges = ranges[:j] + ranges[j + 1:]
                    else:
                        ranges = ranges[:j]
                else:
                    j += 1
            i += 1

        return ranges

    def _filter_forbidden_ranges(self, ranges):

        if self.forbidden_interpro_domains is None:
            raise InitError("forbidden ranges not set")

        forbidden = []
        for range_ in ranges:
            if range_.ac in self.forbidden_interpro_domains:
                forbidden.append(range_)

        passed = []
        for range_ in ranges:
            overlapping = list(filter(lambda r: r.overlaps_with(range_), forbidden))
            if len(overlapping) <= 0:
                passed.append(range_)

        return passed

domain_aligner = DomainAligner()
