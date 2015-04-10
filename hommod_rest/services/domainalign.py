#!/usr/bin/python

import logging

from hommod_rest.services.blast import blaster
from hommod_rest.services.secstr import secstr
from hommod_rest.services.align import aligner

from modelutils import (get_aa321, getNalignIdentity, filterGoodHits,
                        parseDSSP,
                        minIdentity,
                        TemplateID, getTemplatePDBIDandChain)

_log = logging.getLogger(__name__)


def _get_midline(alignedseq1, alignedseq2):

    m = ''
    for i in range(len(alignedseq1)):
        if alignedseq1[i] == alignedseq2[i]:
            m += alignedseq1[i]
        else:
            m += ' '
    return m


def alignmentRepr(aligned, order=[]):

    s = ''

    if len(order) <= 0:
        order = aligned.keys()
    k = order[0]
    for i in range(0, len(aligned[k]), 100):
        f = i + 100
        if f > len(aligned[k]):
            f = len(aligned[k])
        for key in order:
            s += aligned[key][i:f] + '\n'
        s += '\n'

    return s


def _map_gaps(gapped, gapless):
    s = ''
    i = 0
    while len(s) < len(gapped):
        if gapped[len(s)].isalpha():
            s += gapless[i]
            i += 1
        else:
            s += gapped[len(s)]
    return s


def _get_relative_span(alignfrom, alignTo):

    i = 0
    while not alignfrom[i].isalpha():
        i += 1
    start = len(alignTo[:i].replace('-', ''))
    i = len(alignfrom)
    while not alignfrom[i-1].isalpha():
        i -= 1
    end = len(alignTo[:i].replace('-', ''))

    return TargetRange(start, end)


class TargetRange(object):

    def __init__(self, start, end):

        if type(start) != int or type(end) != int:
            raise Exception('only integer ranges allowed')

        self.start = min(start, end)
        self.end = max(start, end)

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

        # given: r[0]>r[1]
        return (self.end >= otherRange.start and self.start <= otherRange.end
                or otherRange.end >= self.start
                and otherRange.start <= self.end)

    def encloses(self, child):

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


class Picker(object):

    def accepts(self, hitID, alignment):

        return True


def distance(range1, range2):  # in AA

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


def overlapStats(range1, range2):

    len1 = range1.length()
    len2 = range2.length()
    noverlap = min(range1.end, range2.end) - max(range1.start, range2.start)
    poverlap = (100.0 * noverlap) / min(len1, len2)
    plendiff = 100.0 * abs(len1 - len2) / max(len1, len2)

    return poverlap, plendiff


def mergeSimilar(ranges):

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

            poverlap, plendiff = overlapStats(ranges[i], ranges[j])

            merged = merge(ranges[i], ranges[j])

            if poverlap > 80.0 and plendiff < 10.0:

                ranges = (ranges[:i] + [merged] + ranges[i + 1:j] +
                          ranges[j + 1:])

        i += 1

        ranges = removeDoubles(ranges)

    ranges.sort()

    return ranges


def removeFromList(subj, toremove):

    r = []
    for x in subj:
        if x not in toremove:
            r.append(x)
    return r


def removeDoubles(ranges):

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


def sortSmallestFirst(ranges):

    sized = []
    for r in ranges:
        sized.append([r.length(), r])
    sized.sort()
    sorted = []
    for l, r in sized:
        sorted.append(r)
    return sorted


def sortLargestFirst(ranges):
    sorted = sortSmallestFirst(ranges)
    sorted.reverse()
    return sorted


def findSharedHits(ranges):

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


class _AlignmentPool(object):  # stores alignments for easy access

    def __init__(self, targetSequence):
        self.pool = {}
        self.targetSequence = targetSequence

    def getAlignment(self, _range, template):

        templateID = str(template)

        if _range not in self.pool:
            self.pool[_range] = {}
        if templateID not in self.pool[_range]:
            self.pool[_range][templateID] = \
                makeDSSPAlignment(
                    self.targetSequence[_range.start: _range.end], template)

        return self.pool[_range][templateID]

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


def hasTemplateSeq(alignment, templateSeq):

    return alignment['template'].replace('-', '') == templateSeq


def makeYasaraAlignment(yasara, domainSeq, yasaraObject, yasaraChainID):

    cas = []
    for s in yasara.ListAtom('CA and protein and obj %i and mol %s' %
                             (yasaraObject, yasaraChainID)):
        cas.append(int(s))

    pdbSeq = ''
    pdbSecStr = ''
    for ca in cas:
        pdbSeq += get_aa321(yasara.ListRes('atom %i' % ca, 'RESNAME')[0])
        pdbSecStr += yasara.SecStrRes('atom %i' % ca)[0]

    aligned = aligner.msaAlign(pdbSeq, pdbSecStr, domainSeq)
    aligned['secstr'] = _map_gaps(aligned['template'], pdbSecStr)
    aligned['midline'] = _get_midline(aligned['template'], aligned['target'])

    if not hasTemplateSeq(aligned, pdbSeq):
        raise Exception(
            'Yasara alignment: aligned seq doesn\'t match template seq:\n' +
            aligned['template'] + '\n' + pdbSeq
        )

    return aligned


def makeDSSPAlignment(domainSeq, template):

    dssppath = '/mnt/cmbi4/dssp/%s.dssp' % template.pdbac.lower()

    d = parseDSSP(dssppath)
    if template.chainID not in d:
        raise Exception('No chain %s in %s' % (template.chainID, dssppath))
    pdbSeq, pdbSecStr, pdbDisulfid = d[template.chainID]

    aligned = aligner.msaAlign(pdbSeq, pdbSecStr, domainSeq)
    aligned['secstr'] = _map_gaps(aligned['template'], pdbSecStr)
    aligned['midline'] = _get_midline(aligned['template'], aligned['target'])

    if not hasTemplateSeq(aligned, pdbSeq):
        raise Exception(
            'Alignment from ' + dssppath +
            ': aligned seq doesn\'t match template seq:\n' +
            aligned['template'] + '\n' + pdbSeq
        )

    return aligned


def getCoveredTargetRange(alignment):

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


def getTargetCoveredRange(alignment, templateseq):

    if templateseq != alignment['template'].replace('-', ''):
        raise Exception('mismatch between alignment and template sequence')

    coveredRangeStart = 0
    while not alignment['target'][coveredRangeStart].isalpha():
        coveredRangeStart += 1
    coveredRangeEnd = coveredRangeStart
    while alignment['target'][coveredRangeEnd].isalpha():
        coveredRangeEnd += 1

    coveredRangeStart = \
        len(alignment['template'][:coveredRangeStart].replace('-', ''))
    coveredRangeEnd = len(templateseq) - \
        len(alignment['template'][coveredRangeEnd:].replace('-', ''))

    return [coveredRangeStart, coveredRangeEnd]


def pickAlignments(yasaraChain, targetSeqs, targetsInterproRanges, picker):

    alignmentTriples = {}
    for targetID in targetsInterproRanges.keys():

        triples = []
        for _range, templateID, alignment in \
                getAlignments(
                    targetsInterproRanges[targetID],
                    targetSeqs[targetID], yasaraChain
                ):

            if 'X' in yasaraChain.seq:
                print alignment['template']

            if picker.accepts(targetID, alignment):

                triples.append((_range, templateID, alignment))

        if len(triples) > 0:
            alignmentTriples[targetID] = triples

    return alignmentTriples


def getForbiddenRanges(interproDomains):

    rs = []
    for domain in interproDomains:
        if domain.ac == 'IPR003596':  # Ig variable domain like
            rs.append(TargetRange(domain.start, domain.end))

    return rs


def getAlignments(interproDomains, tarSeq, yasaraChain=None):

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
        mergedSampleRanges = mergeSimilar(sampleRanges)

        for r in sortLargestFirst(mergedSampleRanges):

            if r in checkedRanges:
                continue
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

            bestTemplate = None
            bestAlignment = None
            bestPID = 0.0
            bestPCOVER = 0.0

            _log.debug('trying range: %s' % r)

            if yasaraChain:

                templateSelected = \
                    TemplateID(yasaraChain.objname, yasaraChain.chainID)

                # align only against 1 template:
                aligned = alignmentDAO.getAlignmentFromYasara(
                    yasaraChain.yasaramodule, r,
                    yasaraChain.obj, yasaraChain.chainID
                )
                nalign, pid = getNalignIdentity(
                    aligned['target'], aligned['template'])
                pcover = (nalign * 100.0) / (r.end - r.start)

                # bHH: highly homologous to the whole target sequence
                bHH = pid >= 80.0 and r.length() == len(tarSeq)

#               print r
#               print aligned['target']
#               print aligned['template']
#               print 'pid:',pid,'minimum:',minIdentity(nalign),'pcover:',pcover
                if pid >= minIdentity(nalign) and pcover >= 80.0 or bHH:

                    if pcover < 80.0:
                        m = getCoveredTargetRange(aligned)
                    else:
                        m = TargetRange(r.start, r.end)
                        checkedRanges.append(m)

                    m.template = templateSelected
                    m.pid = pid
                    m.pcover = pcover
                    okRanges.append(m)

                    if pcover >= 80.0:
                        bestTemplate = templateSelected
                        bestAlignment = aligned
                        bestPID = pid
                        bestPCOVER = pcover

                    _log.debug("domainalign: passing alignment with %s:\n%s"
                               % (templateSelected, alignmentRepr(
                                  aligned, ['target', 'midline', 'template'])))

                else:
                    _log.debug("domainalign: rejecting alignment with %s:\n%s"
                               % (templateSelected, alignmentRepr(
                                  aligned, ['target', 'midline', 'template'])))

            else:

                # iterate over blast hits:
                hits = filterGoodHits(
                    blaster.templateBlast(tarSeq[r.start: r.end]))
                for hitID in hits.keys():

                    pdbid, pdbchain = getTemplatePDBIDandChain(hitID)
                    template = TemplateID(pdbid, pdbchain)
                    if not secstr.hasSecStr(template):

                        _log.warn("domainalign: no secondary structure for %s"
                                  % template)
                        continue

                    _log.debug("domainalign: got blast hit %s" % template)

                    for alignment in hits[hitID]:

                        nalign = alignment.getNumberResiduesAligned()
                        pid = alignment.getIdentity()
                        pcover = (nalign*100.0)/(r.end-r.start)

                        # bHH: highly homologous to the whole target sequence
                        bHH = pid >= 80.0 and r.length() == len(tarSeq)

                        # minimal criteria to make a model:
                        if pid >= minIdentity(nalign) and pcover >= 80.0 or bHH:

                            aligned = {'target': alignment.queryalignment,
                                       'template': alignment.subjectalignment}
                            if pcover < 80.0:
                                m = getCoveredTargetRange(aligned)
                            else:
                                m = TargetRange(r.start, r.end)
                                checkedRanges.append(m)
                            m.template = template
                            m.pid = pid
                            m.pcover = pcover
                            okRanges.append(m)

                            aligned['midline'] = \
                                _get_midline(aligned['target'],
                                             aligned['template'])
                            _log.debug(
                                "domainalign: passing alignment with %s:\n%s"
                                % (template,
                                   alignmentRepr(aligned,
                                                 ['target', 'midline',
                                                  'template'])))

                            # is this one better than previous hits?
                            if pcover >= 80.0 and \
                                    (not bestTemplate or pid > bestPID):

                                bestPID = m.pid
                                bestTemplate = m.template
                                bestAlignment = aligned
                                bestPCOVER = m.pcover

                            break  # we're done with this template

                        else:
                            _log.debug(
                                "domainalign: rejecting blast hit with %s:"
                                % template)

            if bestTemplate:  # we have a best hit for this range

                # Remove any smaller ranges that this one encloses:
                i = 0
                while i < len(bestRanges):
                    if r.encloses(bestRanges[i]):
                        # this one should replace the smaller one
                        bestRanges = bestRanges[: i] + bestRanges[i + 1:]
                    else:
                        i += 1

                aligned = bestAlignment

                if bestPCOVER < 80.0:
                    m = getCoveredTargetRange(aligned)
                else:
                    m = TargetRange(r.start, r.end)
                m.template = bestTemplate
                m.alignment = bestAlignment
                m.pid = bestPID
                m.pcover = bestPCOVER
                bestRanges.append(m)

        # See if we can merge ranges that have
        # the same template in their blast hits:
        checkedRanges = removeDoubles(checkedRanges + sampleRanges)
        sampleRanges = []
        dictSharedHits = findSharedHits(okRanges)
        for template in dictSharedHits.keys():

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

    # populate the returned list:
    returnAlignments = []
    bestRanges.sort()
    for r in bestRanges:

        nalign, pid = getNalignIdentity(aligned['target'], aligned['template'])

        if pid >= minIdentity(nalign):
            returnAlignments.append((r, r.template, r.alignment))

    _log.debug("domainalign: returning %i alignments" % len(returnAlignments))

    return returnAlignments


def alignmentsOfSameTemplateCompatible(alignment1, alignment2):

    span1 = _get_relative_span(alignment1['target'], alignment1['template'])
    span2 = _get_relative_span(alignment2['target'], alignment2['template'])

    # They must not occupy the same piece of template
    return (not span1.overlapsWith(span2))


def swap(x, y):
    return [y, x]


def mergeAlignmentsOfSameTemplate(alignment1, alignment2):

    span1 = _get_relative_span(alignment1['target'], alignment1['template'])
    span2 = _get_relative_span(alignment2['target'], alignment2['template'])

    # make sure span1 is the left sided:
    if span2 < span1:
        span1, span2 = swap(span2, span1)
        alignment1, alignment2 = swap(alignment2, alignment1)

    i1 = 0
    naa = 0
    while naa < span1.end:
        if alignment1['template'][i1].isalpha():
            naa += 1
        i1 += 1

    i2 = 0
    naa = 0
    while naa < span2.start:
        if alignment2['template'][i2].isalpha():
            naa += 1
        i2 += 1

    return {'target': alignment1['target'][: i1] +
            alignment2['target'][: i2],
            'template': alignment1['template'][: i1] +
            alignment2['template'][: i2]}


def joinAlignmentsToBestTemplateCoverage(alignmentTriples):

    # First determine which alignment has the largest cover
    bestAlignment = None
    bestCover = 0.0
    templates = []
    for r, templateID, alignment in alignmentTriples:
        if templateID not in templates:
            if len(templates) > 0:
                raise Exception('not all the same template')
            templates.append(templateID)
        if not bestAlignment or r.pcover > bestCover:
            bestAlignment = alignment.copy()
            bestCover = r.pcover

    # Now merge the best alignment with smaller, compatible alignments
    for r, templateID, alignment in alignmentTriples:
        if alignmentsOfSameTemplateCompatible(bestAlignment, alignment):
            bestAlignment = \
                mergeAlignmentsOfSameTemplate(bestAlignment, alignment)

    return bestAlignment
