#!/usr/bin/python

from modelutils import YasaraChain

import domainalign

import logging

_log = logging.getLogger(__name__)


def getTargetCoveredRange(alignment, templateseq):

    if templateseq != alignment['template'].replace('-', ''):
        raise Exception(
            'mismatch between alignment and template sequence\nseg:' +
            templateseq+'\nali:'+alignment['template']
        )

    coveredRangeStart = 0
    while not alignment['target'][coveredRangeStart].isalpha():
        coveredRangeStart += 1
    coveredRangeEnd = coveredRangeStart
    while (coveredRangeEnd < len(alignment['target']) and
           alignment['target'][coveredRangeEnd].isalpha()):

        coveredRangeEnd += 1

    coveredRangeStart = len(alignment['template'][:coveredRangeStart]
                            .replace('-', ''))
    coveredRangeEnd = (len(templateseq) -
                       len(alignment['template'][coveredRangeEnd:]
                           .replace('-', '')))

    return [coveredRangeStart, coveredRangeEnd]


def listInteractingChains(yasaraChain):

    return yasaraChain.yasaramodule.ListMol(
        'protein and obj %i and not mol %s with distance<4.5 from obj %i mol %s'
        % (yasaraChain.obj, yasaraChain.chainID,
           yasaraChain.obj, yasaraChain.chainID), 'MOL')

    chains = []

    for ca in yasaraChain.CAs:
        for chain2 in (yasaraChain.yasaramodule.ListAtom(
                       'CA and obj %i with distance<6 from %i' %
                       (yasaraChain.obj, ca), 'MOL')):
            if chain2 != yasaraChain.chainID and chain2 not in chains:
                chains.append(chain2)

    return chains


class InteractionPicker(domainalign.Picker):

    def __init__(self, subjectYasaraChain, interactionChainAlignments):
        # Alignments tell which regions are covered

        self.subjectYasaraChain = subjectYasaraChain
        self.interactionChainAlignments = interactionChainAlignments

    def accepts(self, targetID, subjectChainAlignment):

        for chainID in self.interactionChainAlignments.keys():

            sCAs = self.subjectYasaraChain.CAs
            sseq = self.subjectYasaraChain.seq

            interactionYasaraChain = YasaraChain(
                self.subjectYasaraChain.yasaramodule,
                self.subjectYasaraChain.obj, chainID
            )
            iCAs = interactionYasaraChain.CAs
            iseq = interactionYasaraChain.seq

            sstart, send = getTargetCoveredRange(subjectChainAlignment, sseq)
            istart, iend = getTargetCoveredRange(
                self.interactionChainAlignments[chainID], iseq
            )

            for i in range(sstart, send):
                sca = sCAs[i]

                search_cmd = ('CA and atom {}-{} and obj {} and ' +
                              'mol {} with distance<6 from {}') \
                    .format(iCAs[istart], iCAs[iend - 1],
                            interactionYasaraChain.obj,
                            interactionYasaraChain.chainID, sca)

                interactingCAs =\
                    interactionYasaraChain.yasaramodule.ListAtom(search_cmd)

                if 0 < len(interactingCAs):

                    _log.debug("found interaction between chains {} and {}"
                               .format(interactionYasaraChain.chainID,
                                       self.subjectYasaraChain.chainID))

                    return True
            # Two atoms, covered in the two alignments, are close to each other

        return False
