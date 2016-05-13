#!/usr/bin/python

# The code in this script was added to account for interactions
# between protein chains in the template. When making an alignment,
# We'd like to preserve as much interaction as possible...

from modelutils import YasaraChain

import domainalign

import logging

_log = logging.getLogger(__name__)


# This functino must return the start and end position of the range
# in the target sequence, covered by the given template.
# Assumes that the alignment has the sequence ids 'target' and 'template'
# Alignment target sequence is assumed to be the full target sequence.
def getTargetCoveredRange (alignment, templateseq):

    # templateseq must match its alignment ['template'],
    # but alignment ['template'] may be a substring of templateseq:
    if templateseq != alignment ['template'].replace('-', ''):
        raise Exception(
            'mismatch between alignment and template sequence\nseq:' +
            templateseq+'\nali:'+alignment['template']
        )

    # Determine alignment['target']'s start and end positions in alignment:
    targetStart = 0
    while not alignment['target'][targetStart].isalpha():

        targetStart += 1

    targetEnd = targetStart
    while (targetEnd < len(alignment['target']) and
           alignment['target'][targetEnd].isalpha()):

        targetEnd += 1


    # Determine where alignment['template'] starts, relative to alignment['target']:
    coveredRangeStart = len(alignment['template'][:targetStart]
                            .replace('-', ''))

    coveredRangeEnd = (len(templateseq) -
                       len(alignment['template'][targetEnd:]
                           .replace('-', '')))

    return (coveredRangeStart, coveredRangeEnd)


# This function must tell which tell which other template chains
# interact with the given one.
def listInteractingChains (yasaraChain):

    # Take everything that's close enough:
    return yasaraChain.yasaramodule.ListMol(
        'protein and obj %i and not mol %s with distance<4.5 from obj %i mol %s'
        % (yasaraChain.obj, yasaraChain.chainID,
           yasaraChain.obj, yasaraChain.chainID), 'MOL')

#    chains = []
#
#    for ca in yasaraChain.CAs:
#        for chain2 in (yasaraChain.yasaramodule.ListAtom(
#                       'CA and obj %i with distance<6 from %i' %
#                       (yasaraChain.obj, ca), 'MOL')):
#            if chain2 != yasaraChain.chainID and chain2 not in chains:
#                chains.append(chain2)
#
#    return chains


# This class accepts/rejects alignments, based on interaction preservation.
# On the basis of the subject chain and other chains that interact with it.
# Takes alignments for both chains and checks whether the target-covered parts
# are interacting in those alignments.
class InteractionPicker (domainalign.Picker):

    def __init__ (self, subjectChainID, yasaraChains, interactionChainAlignments):
        # Alignments tell which regions are covered

        self.yasaraChains = yasaraChains
        self.subjectChainID = subjectChainID

        if subjectChainID not in yasaraChains:
            raise Exception ("must provide all yasara chains of the complex, " +
                             "including subject chain %s" % subjectChainID)

        self.interactionChainAlignments = interactionChainAlignments

        for chainID in interactionChainAlignments:
            if chainID != subjectChainID:
                if chainID not in yasaraChains:
                    raise Exception ("yasara chain not provided for %s"
                                     % chainID)

    def accepts (self, targetID, subjectChainAlignment):

        # Check all interacting chains' alignments for interactions:
        for chainID in self.interactionChainAlignments:

            # Determine which residues in the chains are covered by a target sequence:
            subjectYasaraChain = self.yasaraChains [self.subjectChainID]
            sCAs = subjectYasaraChain.CAs
            sseq = subjectYasaraChain.seq

            interactionYasaraChain = self.yasaraChains [chainID]

            iCAs = interactionYasaraChain.CAs
            iseq = interactionYasaraChain.seq

            sstart, send = getTargetCoveredRange(subjectChainAlignment, sseq)
            istart, iend = getTargetCoveredRange(
                self.interactionChainAlignments[chainID], iseq
            )

            # Check every target-covered residue in the subject chain:
            for i in range (sstart, send):

                # Subject chain C-alpha:
                sca = sCAs [i]

                # Search for close C-alphas in the other chain:
                search_cmd = ('CA and atom {}-{} and obj {} and ' +
                              'mol {} with distance<6 from {}') \
                    .format(iCAs[istart], iCAs[iend - 1],
                            interactionYasaraChain.obj,
                            interactionYasaraChain.chainID, sca)

                interactingCAs =\
                    interactionYasaraChain.yasaramodule.ListAtom(search_cmd)

                # Return True if a single interacting C-alpha pair is found:
                if 0 < len (interactingCAs):

                    _log.debug("found interaction between chains {} and {}"
                               .format(interactionYasaraChain.chainID,
                                       subjectYasaraChain.chainID))

                    return True
            # Two atoms, covered in the two alignments, are close to each other

        return False
