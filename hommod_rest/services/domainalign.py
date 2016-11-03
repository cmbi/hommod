#!/usr/bin/python

import logging
import traceback

from hommod_rest.services.blast import blaster
from hommod_rest.services.secstr import secstr
from hommod_rest.services.align import aligner

from modelutils import (get_aa321, getNalignIdentity, filterGoodHits,
                        parseDSSP,
                        minIdentity,
                        TemplateID, getTemplatePDBIDandChain)

_log = logging.getLogger(__name__)

SIMILAR_RANGES_MIN_OVERLAP_PERCENT = 80.0
SIMILAR_RANGES_MAX_LENDIFF_PERCENT = 10.0

# Loads the template blacklist file and returns whether the given
# pdb id occurs in it or not.
# (blacklist file is set in configuration)
def _template_in_blacklist (pdbid):

    from hommod_rest.services.model import modeler
    return modeler._template_in_blacklist (pdbid)

# This generates the midline string for an alignment.
# It's used for debugging purposes.
def _get_midline (alignedseq1, alignedseq2):

    m = ''
    for i in range (len (alignedseq1)):
        if alignedseq1 [i] == alignedseq2 [i]:
            m += alignedseq1 [i]
        else:
            m += ' '
    return m

# Represents a given alignment dictionary as a string.
# This is used for debugging.
def alignmentRepr (aligned, order=[]):

    s = ''

    if len(order) <= 0:
        order = aligned.keys()
    k = order[0]
    for i in range(0, len(aligned[k]), 100):
        f = i + 100
        if f > len (aligned[k]):
            f = len (aligned[k])
        for key in order:
            s += aligned[key][i:f] + '\n'
        s += '\n'

    return s


# maps the sequence of 'gapless' onto
# 'gapped', keeping the gaps, but replacing the
# amino acids in sequential order.
def _map_gaps (gapped, gapless):

    s = ''
    i = 0
    while len (s) < len (gapped):

        if gapped [len (s)].isalpha():
            s += gapless [i]
            i += 1
        else:
            s += gapped [len (s)]
    return s


# Assumes that 'alignfrom' and 'alignto' are both
# sequences (strings) in the same alignment.
# It tells the starting position of 'alignfrom'
# relative to 'the starting position of 'alignTo'
def _get_relative_span (alignfrom, alignTo):

    i = 0
    while not alignfrom[i].isalpha():
        i += 1
    start = len(alignTo[:i].replace('-', ''))
    i = len(alignfrom)
    while not alignfrom[i-1].isalpha():
        i -= 1
    end = len(alignTo[:i].replace('-', ''))

    return TargetRange (start, end)


# The following class represents a range on a target
# sequence. Most important fields are the start and end
# positions. These are assumed to be ranging from 0 to
# the length of the sequence.
class TargetRange(object):

    def __init__(self, start, end):

        if type(start) != int or type(end) != int:
            raise Exception('only integer ranges allowed')

        self.start = min(start, end)
        self.end = max(start, end)

    # If range 1 lies left of range 2, then range1 < range2
    def __lt__(self, other):

        if self.start < other.start:
            return True
        elif self.start == other.start:
            return self.end < other.end
        else:
            return False

    def __eq__(self, other):

        return (self.start == other.start and self.end == other.end)

    def __hash__(self):

        return hash((self.start, self.end))

    def __repr__(self):

        return "%i-%i" % (self.start, self.end)

    def length(self):

        return (self.end - self.start)

    def overlapsWith(self, otherRange):

        # given: start <= end
        return (self.end >= otherRange.start and self.start <= otherRange.end
                or otherRange.end >= self.start
                and otherRange.start <= self.end)

    def encloses(self, child): # all of 'child' lies in 'self'

        return (self.start <= child.start and self.end >= child.end)

    def findEnclosedChildren(self, candidateRanges):
        indices = []

        for i in range(len(candidateRanges)):
            if self.encloses(candidateRanges[i]):
                indices.append(i)

        return indices

    def findOverlapping(self, candidateRanges):
        indices = []

        for i in range(len(candidateRanges)):
            if self.overlapsWith(candidateRanges[i]):
                indices.append(i)

        return indices

    def __add__ (self, rshift):

      return TargetRange (self.start + rshift, self.end + rshift)

# A picker object can filter alignments, when inserted
# into the pickAlignment function.
# (must be implemented elsewhere)
class Picker (object):

    def accepts (self, hitID, alignment):

        return True


def distance (range1, range2):  # in AA

    if range1.start < range2.start:

        if range1.end < range2.start:

            return range2.start - range1.end

        else:  # overlap

            return 0

    else:  # range2.start < range1.start

        if range2.end < range1.start:

            return range1.start - range2.end

        else:  # overlap

            return 0


def merge(range1, range2):

    return TargetRange(min(range1.start, range2.start),
                       max(range1.end, range2.end))


def intersection(range1, range2):

    return TargetRange(max(range1.start, range2.start),
                       min(range1.end, range2.end))

# Returns the % overlap between the two ranges and the
# % of length difference, relative to the largest of the two ranges.
def overlapStats(range1, range2):

    len1 = range1.length()
    len2 = range2.length()
    noverlap = min(range1.end, range2.end) - max(range1.start, range2.start)
    poverlap = (100.0 * noverlap) / min(len1, len2)
    plendiff = 100.0 * abs(len1 - len2) / max(len1, len2)

    return poverlap, plendiff

# Merges ranges of significant overlap, with almost the same length.
# Returns the list of merged ranges.
def mergeSimilar (ranges):

    ranges.sort()

    # Join everything that's highly similar:
    i = 0
    while i < len(ranges):
        overlapping = ranges[i].findOverlapping(ranges)
        overlapping.sort()
        # important, largest index must go first!
        # Because we're goin to remove ranges from the list.
        overlapping.reverse()

        for j in overlapping:

            if j == i:
                continue  # don't merge with itself

            poverlap, plendiff = overlapStats (ranges[i], ranges[j])

            merged = merge(ranges[i], ranges[j])

            if poverlap > SIMILAR_RANGES_MIN_OVERLAP_PERCENT and \
                    plendiff < SIMILAR_RANGES_MAX_LENDIFF_PERCENT:

                ranges = (ranges[:i] + [merged] + ranges[i + 1:j] +
                          ranges[j + 1:])

        i += 1

        # Make list shorter to save time:
        ranges = removeDoubles (ranges)

    ranges.sort()

    return ranges

# Returns the given list without the elements in 'toremove'.
def removeFromList(subj, toremove):

    r = []
    for x in subj:
        if x not in toremove:
            r.append(x)
    return r

# Returns the given list, with doubles removed.
def removeDoubles (ranges):

    i = 0
    while i < len(ranges):
        j = i + 1
        while j < len(ranges):
            if ranges[i] == ranges[j]:

                if (j + 1) < len (ranges):
                    ranges = ranges[:j] + ranges[j + 1:]
                else:
                    ranges = ranges[:j]
            else:
                j += 1
        i += 1

    return ranges

# Returns the given list with smallest ranges in front.
def sortSmallestFirst (ranges):

    sized = []
    for r in ranges:
        sized.append([r.length(), r])
    sized.sort()
    sorted = []
    for l, r in sized:
        sorted.append(r)
    return sorted

# Returns the given list with largest ranges in front.
def sortLargestFirst (ranges):

    sorted = sortSmallestFirst(ranges)
    sorted.reverse()
    return sorted

# Looks for ranges that have the 'template' field set.
# If two ranges have the same template, they are put
# together in the returned dictionary as a list.
def findSharedHits (ranges):

    # sort by template id
    d = {}
    for r in ranges:
        if hasattr(r, 'template'):
            if r.template not in d:
                d[r.template] = []
            d[r.template].append(r)

    # Forget those blast hits that only occur once:
    r = {}
    for template in d.keys():
        if len(d[template]) > 1:
            r[template] = d[template]

    return r


# stores alignments for easy access.
# Lookup alignment by range object and template (structure + chain id)
class _AlignmentPool(object):

    def __init__(self, targetSequence):
        self.pool = {}
        self.targetSequence = targetSequence

    # For lookup by pdbid:
    def getAlignment(self, _range, template):

        templateID = str(template)

        if _range not in self.pool:
            self.pool[_range] = {}
        if templateID not in self.pool[_range]:
            self.pool[_range][templateID] = \
                makeDSSPAlignment(
                    self.targetSequence[_range.start: _range.end], template)

        return self.pool[_range][templateID]

    # For lookup by yasara object:
    def getAlignmentFromYasara(self, yasara, _range, yasaraObj, yasaraChain):

        templateID = 'y-' + \
                     str(TemplateID(
                         yasara.ListObj(str(yasaraObj),
                                        'OBJNAME')[0], yasaraChain))

        if _range not in self.pool:
            self.pool[_range] = {}
        if templateID not in self.pool[_range]:
            self.pool[_range][templateID] = \
                makeYasaraAlignment(
                    yasara,
                    self.targetSequence[_range.start:_range.end],
                    yasaraObj, yasaraChain
                )

        return self.pool[_range][templateID]

# Returns true if the given templateSeq matches the alignment's
# template sequence, without gaps.
def hasTemplateSeq (alignment, templateSeq):

    return alignment['template'].replace('-', '') == templateSeq


# Produces an alignment from a yasara object as template seq.
def makeYasaraAlignment (yasara, domainSeq, yasaraObject, yasaraChainID):

    cas = []
    for s in yasara.ListAtom('CA and protein and obj %i and mol %s' %
                             (yasaraObject, yasaraChainID)):
        cas.append(int(s))

    pdbSeq = ''
    pdbSecStr = ''
    for ca in cas:
        pdbSeq += get_aa321(yasara.ListRes('atom %i' % ca, 'RESNAME')[0])
        pdbSecStr += yasara.SecStrRes('atom %i' % ca)[0]

    aligned = aligner.kmad_align(pdbSeq, pdbSecStr, domainSeq)
    aligned['secstr'] = _map_gaps(aligned['template'], pdbSecStr)
    aligned['midline'] = _get_midline(aligned['template'], aligned['target'])

    if not hasTemplateSeq(aligned, pdbSeq):
        raise Exception(
            'Yasara alignment: aligned seq doesn\'t match template seq:\n' +
            aligned['template'] + '\n' + pdbSeq
        )

    return aligned

# Produces an alignment, using the dssp file as template.
def makeDSSPAlignment(domainSeq, template):

    dssppath = '/mnt/cmbi4/dssp/%s.dssp' % template.pdbac.lower()

    d = parseDSSP(dssppath)
    if template.chainID not in d:
        raise Exception('No chain %s in %s' % (template.chainID, dssppath))
    pdbSeq, pdbSecStr, pdbDisulfid = d[template.chainID]

    aligned = aligner.kmad_align(pdbSeq, pdbSecStr, domainSeq)
    aligned['secstr'] = _map_gaps(aligned['template'], pdbSecStr)
    aligned['midline'] = _get_midline(aligned['template'], aligned['target'])

    if not hasTemplateSeq(aligned, pdbSeq):
        raise Exception(
            'Alignment from ' + dssppath +
            ': aligned seq doesn\'t match template seq:\n' +
            aligned['template'] + '\n' + pdbSeq
        )

    return aligned

# Give this function an alignment with a 'target' and 'template'
# sequence, then it will return the range, representing the
# part of the target, that is covered by the template.
def getCoveredTargetRange (alignment):

    targetStart = 0
    while not alignment['target'][targetStart].isalpha():
        targetStart += 1

    targetEnd = len(alignment['target']) - 1
    while not alignment['target'][targetEnd].isalpha():
        targetEnd -= 1

    start = targetStart
    aaStart = 0
    while not alignment['template'][start].isalpha():
        if alignment['target'][start].isalpha():
            aaStart += 1
        start += 1

    end = targetEnd
    aaEnd = 0
    while not alignment['template'][end].isalpha():
        if alignment['target'][end].isalpha():
            aaEnd += 1
        end -= 1

    return TargetRange(targetStart + aaStart, targetEnd - aaEnd)

# For the given target sequence range, this function
# returns the corresponding sequence of the template,
# according to the given alignment.
def getTemplateSeqAtTargetPositions(alignment, startInTarget, endInTarget):

    targetStart = 0
    while not alignment['target'][targetStart].isalpha():
        targetStart += 1

    start = targetStart
    nAA = 0
    while nAA < startInTarget and start < len(alignment['target']):
        if alignment['target'][start].isalpha():
            nAA += 1
        start += 1

    end = start
    nAA = 0
    while nAA < (endInTarget - startInTarget) and \
            end < len(alignment['target']):

        if alignment['target'][end].isalpha():
            nAA += 1
        end += 1

    return alignment['template'][start: end]

# For a given alignment and matching template sequence, this function
# returns the corresponding template range that is covered by the target.
def getTargetCoveredRange (alignment, templateseq):

    if templateseq != alignment['template'].replace('-', ''):
        raise Exception('mismatch between alignment and template sequence')

    coveredRangeStart = 0
    while not alignment['target'][coveredRangeStart].isalpha():
        coveredRangeStart += 1
    coveredRangeEnd = coveredRangeStart
    while alignment['target'][coveredRangeEnd].isalpha():
        coveredRangeEnd += 1

    coveredRangeStart = \
        len (alignment['template'][:coveredRangeStart].replace('-', ''))
    coveredRangeEnd = len(templateseq) - \
        len (alignment['template'][coveredRangeEnd:].replace('-', ''))

    return (coveredRangeStart, coveredRangeEnd)


# Generates alignments using the 'getAlignments' function
# and filters out those, approved by the inserted picker object.
# It runs the 'getAlignments' function on the given targets and template.
# It's typically used when you already have a template complex, but need to find
# targets for some template chains, other than the main target's template chain.
#
# The returned dictionary of triples uses the same keys as the inserted
# 'targetSeqs' and 'targetsInterproRanges' dictionaries. (must match)
# The returned dictionary's values are triples containing the template covered
# range, the template ID (pdbid and chain) and the target-template alignment.
def pickAlignments (templateYasaraChain, targetSeqs, targetsInterproRanges, picker):

    alignmentTriples = {}
    for targetID in targetsInterproRanges.keys():

        triples = []
        for _range, templateID, alignment in \
                getAlignments(
                    targetsInterproRanges[targetID],
                    targetSeqs[targetID], templateYasaraChain
                ):

            if picker.accepts (targetID, alignment):

                triples.append((_range, templateID, alignment))

        if len(triples) > 0:
            alignmentTriples[targetID] = triples

    return alignmentTriples

# We don't want to model Immunoglobulin variable chains.
# So This function must point out ranges that represent those.
# For now, we only have one interpro entry blacklisted.
def getForbiddenRanges (interproDomains):

    rs = []
    for domain in interproDomains:
        if domain.ac == 'IPR003596':  # Ig variable domain like
            rs.append(TargetRange(domain.start, domain.end))

    return rs

def _alignment_ok_for_range (r, tarSeq, nalign, pid, pcover):

    # bHH: highly homologous to the whole target sequence
    bHH = pid >= 80.0 and r.length() == len (tarSeq)

    return pid >= minIdentity (nalign) and (pcover >= 80.0 or bHH)

def _get_range_from (r, template, alignment, pid, pcover):

    if pcover < 80.0:
        m = getCoveredTargetRange(alignment)
    else:
        m = TargetRange(r.start, r.end)

    m.template = template
    m.pid = pid
    m.pcover = pcover

    return m

class _Best_Hit_candidate (object):

    def __init__ (self, template, alignment, pid, pcover):

        self.template = template
        self.alignment = alignment
        self.pid = pid
        self.pcover = pcover

# This function tries to find the best possible sets of alignments,
# given the target sequence, interpro ranges and optionally, a template object.
# * Joins overlapping ranges
# * Removes ranges that are similar
# * Throws out alignments, completely enclosed by a larger one.
# * Checks for sufficient coverage and identity.
# * If a yasara chain is given, it will only make alignments with that template.
def getAlignments (interproDomains, tarSeq, yasaraChain=None):

    _log.debug("get alignments on %s" % tarSeq)

    forbiddenRanges = getForbiddenRanges(interproDomains)

    # Get the domain ranges from interpro:
    # (skip the ones that hit forbidden ranges)
    sampleRanges = []
    for domain in interproDomains:

        r = TargetRange(domain.start, domain.end)
        if len(r.findOverlapping(forbiddenRanges)) <= 0:

            sampleRanges.append(TargetRange(domain.start, domain.end))

    # Add the whole sequence too,
    # just in case the interpro file has no domains for us
    sampleRanges.append(TargetRange(0, len(tarSeq)))

    okRanges = []  # ranges that passed the criteria
    bestRanges = []  # ranges chosen to be modeled

    alignmentDAO = _AlignmentPool(tarSeq)  # allows to access alignments

    checkedRanges = []  # keeps track of all the ranges that we've passed by
    while len(sampleRanges) > 0:

        _log.debug("iter with %i sample ranges" % len(sampleRanges))

        # merging must reduce the number of blasts/alignments:
        mergedSampleRanges = mergeSimilar (sampleRanges)

        # Check the largest ranges first. If that yields, then the smaller
        # ones don't matter.
        for r in sortLargestFirst (mergedSampleRanges):

            if r in checkedRanges:
                continue # Already passed this one

            checkedRanges.append(r)

            # Check to see if this range is part of a larger range that also
            # passed the criteria, if so then don't bother blasting this one
            hasLarger = False
            for i in range(len(bestRanges)):
                if bestRanges[i].encloses(r):
                    hasLarger = True
                    break
            if hasLarger:
                continue

            # Template to pick when there are good choices
            bestHitForRange = None

            # Template to pick when there's only a highly homologous one, with
            # low coverage
            lastResortForRange = None

            _log.debug ('trying range: %s' % r)

            if yasaraChain: # template already chosen outside this function

                templateSelected = \
                        TemplateID(yasaraChain.objname[:4], yasaraChain.chainID)

                # align only against 1 template:
                try:
                    aligned = alignmentDAO.getAlignmentFromYasara(
                        yasaraChain.yasaramodule, r,
                        yasaraChain.obj, yasaraChain.chainID
                    )
                except Exception as e:
                    _log.error("alignment failed:\n" + traceback.format_exc())
                    # If kmad fails, then skip this one :(
                    continue

                nalign, pid = getNalignIdentity(
                    aligned['target'], aligned['template'])
                pcover = (nalign * 100.0) / (r.end - r.start)

                if _alignment_ok_for_range (r, tarSeq, nalign, pid, pcover):

                    # This range made an OK alignment, so at least store it for
                    # later usage:
                    m = _get_range_from (r, templateSelected, aligned,
                                     pid, pcover)
                    okRanges.append (m)

                    # Remember to use this alignment for this range:
                    bestHitForRange = _Best_Hit_candidate (templateSelected,
                                                           aligned,
                                                           pid, pcover)

                    _log.debug ("domainalign: passing %d - %d alignment with %s:\n%s"
                                % (r.start, r.end, templateSelected, alignmentRepr(
                                    aligned, ['target', 'midline', 'template'])))

                else:
                    _log.debug (("domainalign: rejecting %d - %d alignment with %s:\n"
                                 + "%s\ncover: %.1f %%, identity: %.1f %%\n")
                                 % (r.start, r.end, templateSelected, alignmentRepr (
                                        aligned, ['target', 'midline', 'template']).strip (), pcover, pid))

            else: # No chosen template, check all

                _log.debug ('blasting range: %s\n%s' %
                                (r, tarSeq [r.start:r.end]))

                # iterate over blast hits:
                hits = filterGoodHits(
                    blaster.templateBlast(tarSeq[r.start: r.end]))

                _log.debug ('iterating %d blast hits' % len (hits))
                for hitID in hits:

                    pdbid, pdbchain = getTemplatePDBIDandChain(hitID)
                    if _template_in_blacklist (pdbid):

                        _log.debug ("not using %s as template" % pdbid +
                                    ", because it's blacklisted")
                        continue

                    template = TemplateID(pdbid, pdbchain)
                    if not secstr.hasSecStr(template):

                        _log.warn("domainalign: no secondary structure for %s"
                                  % template)
                        continue

                    _log.debug("domainalign: got blast hit %s" % template)

                    for alignment in hits[hitID]:

                        aligned = {'target': alignment.queryalignment,
                                   'template': alignment.subjectalignment}
                        aligned ['midline'] = _get_midline (aligned ['target'],
                                                            aligned ['template'])

                        nalign = alignment.getNumberResiduesAligned ()
                        pid = alignment.getIdentity ()
                        pcover = (nalign * 100.0) / (r.end - r.start)

                        # minimal criteria to make a model:
                        if _alignment_ok_for_range (r, tarSeq, nalign, pid, pcover):

                            # This range made an OK alignment, so at least
                            # store it for later usage:
                            m = _get_range_from (r, template,
                                                 aligned, pid, pcover)
                            okRanges.append(m)

                            _log.debug (
                                "domainalign: passing %d - %d alignment with %s:\n%s"
                                % (r.start, r.end, template,
                                   alignmentRepr (aligned,
                                                  ['target', 'midline',
                                                   'template'])))

                            # is this one better than previous hits?
                            if pcover >= 80.0 and \
                                    (not bestHitForRange or 
                                     pid > bestHitForRange.pid):

                                # Remember to use this alignment for this
                                # range:
                                bestHitForRange = _Best_Hit_candidate (template,
                                                                       aligned,
                                                                       m.pid,
                                                                       m.pcover)
                            else:
                                # Only use low-coverage alignments if no other
                                # options.
                                lastResortForRange = _Best_Hit_candidate (m.template,
                                                                          aligned,
                                                                          m.pid,
                                                                          m.pcover)
                            # we're done with this hitID, since we made
                            # a full alignment.
                            break
                        else:
                            _log.debug(
                                ("domainalign: rejecting %d - %d blast hit with %s:\n"
                                 + "%s\ncover: %.1f %%, identity: %.1f %%\n")
                                 % (r.start, r.end, template,
                                    alignmentRepr (aligned,
                                                   ['target', 'midline',
                                                    'template']).strip (), pcover, pid))


            if not bestHitForRange and lastResortForRange:
                bestHitForRange = lastResortForRange

            if bestHitForRange:  # we have a best hit for this range

                _log.debug ("pick %d - %d best hit %s:\n%s"
                            % (r.start, r.end, bestHitForRange.template,
                               alignmentRepr (bestHitForRange.alignment,
                                              ['target', 'midline',
                                               'template'])))

                # Remove any smaller ranges that this one encloses:
                bestRanges = _remove_enclosing (r, bestRanges)

                # Memorize this best hit:
                m = _range_of_alignment_in (bestHitForRange.alignment, tarSeq)
                m.template = bestHitForRange.template
                m.alignment = bestHitForRange.alignment
                m.pid = bestHitForRange.pid
                m.pcover = bestHitForRange.pcover
                bestRanges.append(m)

        # After iterating the sample ranges once, prepare for the next round:
        sampleRanges = _clean_ranges_search_space (checkedRanges,
                                                   sampleRanges,
                                                   okRanges)


    # Removed the metaphoric question marks from all the sample ranges.
    # Populate the returned list of alignments to model with:
    returnAlignments = []
    bestRanges.sort()
    for r in bestRanges:

        nalign, pid = getNalignIdentity(aligned['target'], aligned['template'])

        if pid >= minIdentity (nalign):

            returnAlignments.append((r, r.template, r.alignment))
        else:
            _log.debug ("throwing away earlier matched %d-%d with %s, because of low identity" %
                        (r.start, r.end, r.template))

    _log.debug ("domainalign: returning %i alignments" % len(returnAlignments))

    return returnAlignments

def _range_of_alignment_in (alignment, tarSeq):

    # We need to translate the alignment's range to the right,
    # because it's only a fraction of the complete target sequence.

    coverstart = tarSeq.find (alignment ['target'].replace ('-',''))
    if coverstart < 0:
        raise Exception ("aligned target sequence not " +
                         "found in full target sequence")

    return getCoveredTargetRange (alignment) + coverstart

def _remove_enclosing (r, ranges):

    i = 0
    while i < len (ranges):
        if r.encloses (ranges [i]):
            _log.debug ("removing smaller %d - %d, " % 
                        (ranges[i].start, ranges [i].end) +
                        " enclosed by larger %d - %d" %
                        (r.start, r.end))

            # this one should replace the smaller one
            ranges = ranges[: i] + ranges[i + 1:]
        else:
            i += 1

    return ranges

def _clean_ranges_search_space (checkedRanges, sampleRanges, okRanges):

    # See if we can merge ranges that have
    # the same template in their blast hits:
    checkedRanges = removeDoubles (checkedRanges + sampleRanges)
    sampleRanges = []
    dictSharedHits = findSharedHits (okRanges)
    for template in dictSharedHits:

        ranges = dictSharedHits[template]

        for i in range(len(ranges)):
            overlapping = ranges[i].findOverlapping(ranges)

            for j in overlapping:
                if j == i:
                    continue  # don't merge with itself

                poverlap, plendiff = overlapStats(ranges[i], ranges[j])
                dist = distance(ranges[i], ranges[j])

                merged = merge(ranges[i], ranges[j])

                # Merge only if:
                # - the ranges are close together
                # - the merge has not already been done
                # - the intersecting parts of the ranges align to
                #   the template in exactly the same way
                if dist < 10 and merged not in checkedRanges:

                    try:
                        if yasaraChain and \
                                yasaraChain.getTemplateID() == template:
                            alignment1 = alignmentDAO.getAlignmentFromYasara(
                                ranges[i], yasaraChain.obj, yasaraChain.chainID)
                            alignment2 = alignmentDAO.getAlignmentFromYasara(
                                ranges[j], yasaraChain.obj, yasaraChain.chainID)
                            alignmentM = alignmentDAO.getAlignmentFromYasara(
                                merged, yasaraChain.obj, yasaraChain.chainID)
                        else:
                            alignment1 = alignmentDAO.getAlignment(
                                ranges[i], template)
                            alignment2 = alignmentDAO.getAlignment(
                                ranges[j], template)
                            alignmentM = alignmentDAO.getAlignment(
                                merged, template)
                    except:
                        # If kmad fails, then skip this one :(
                        continue

                    intersected = intersection(ranges[i], ranges[j])
                    isectTempSeq1 = getTemplateSeqAtTargetPositions(
                        alignment1, intersected.start - ranges[i].start,
                        intersected.end - ranges[i].start)
                    isectTempSeq2 = getTemplateSeqAtTargetPositions(
                        alignment2, intersected.start - ranges[j].start,
                        intersected.end - ranges[j].start)
                    isectTempSeqM = getTemplateSeqAtTargetPositions(
                        alignmentM, intersected.start - merged.start,
                        intersected.end - merged.start)

                    if isectTempSeq1 == isectTempSeqM and \
                            isectTempSeq2 == isectTempSeqM:
                        sampleRanges.append(merged)

    return sampleRanges

# This function checks to be sure tha two alignments are compatible
# i.e. they don't overlap.
def alignmentsOfSameTemplateCompatible (alignment1, alignment2):

    span1 = _get_relative_span (alignment1 ['target'], alignment1 ['template'])
    span2 = _get_relative_span (alignment2 ['target'], alignment2 ['template'])

    # They must not occupy the same piece of template
    return (not span1.overlapsWith (span2))


def swap (x, y):
    return (y, x)


# This function merges two alignments that involve the same
# target and template, but different pieces of them.
def mergeAlignmentsOfSameTemplate (alignment1, alignment2):

    span1 = _get_relative_span (alignment1 ['target'], alignment1 ['template'])
    span2 = _get_relative_span (alignment2 ['target'], alignment2 ['template'])

    # make sure span1 is the left sided:
    if span2 < span1:
        span1, span2 = swap (span1, span2)
        alignment1, alignment2 = swap (alignment1, alignment2)

    # Determine alignment positions of span1 and span2:
    i1 = 0
    naa = 0
    while naa < span1.end:
        if alignment1 ['template'][i1].isalpha():
            naa += 1
        i1 += 1

    i2 = 0
    naa = 0
    while naa < span2.start:
        if alignment2 ['template'][i2].isalpha():
            naa += 1
        i2 += 1

    # Merge alignments:
    return {'target': alignment1['target'][: i1] +
            alignment2 ['target'][: i2],
            'template': alignment1['template'][: i1] +
            alignment2 ['template'][: i2]}

# The following function merges alignments if 'alignmentsOfSameTemplateCompatible'
# returns true for them. It tries to make the merged result as large as possible.
# 'alignmentTriples' is the output from 'getAlignments'
def joinAlignmentsToBestTemplateCoverage (alignmentTriples):

    # First determine which alignment has the largest cover
    bestAlignment = None
    bestCover = 0.0
    templates = []
    for r, templateID, alignment in alignmentTriples:
        if templateID not in templates:
            if len (templates) > 0:
                raise Exception('not all the same template')
            templates.append (templateID)
        if not bestAlignment or r.pcover > bestCover:
            bestAlignment = alignment.copy()
            bestCover = r.pcover

    # Now merge the best alignment with smaller, compatible alignments
    for r, templateID, alignment in alignmentTriples:
        if alignmentsOfSameTemplateCompatible(bestAlignment, alignment):
            bestAlignment = \
                mergeAlignmentsOfSameTemplate(bestAlignment, alignment)

    return bestAlignment
