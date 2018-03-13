import logging
import os
import re
import shutil
import sys
import tarfile
import traceback
from gzip import GzipFile
from StringIO import StringIO

from time import time
from suds import client

from modelutils import (getChainCAsSeqSecStr, getNalignIdentity, downloadPDB,
                        getCoverageIdentity, YasaraChain,
                        identifyDeletedRegions, getPercentageIdentity,
                        filterMinLength, minIdentity, get_pdb_contents,
                        getUniprotSeq, idForSeq, getTemplateSequence)

from interaction import listInteractingChains, InteractionPicker

from filelock import FileLock
from hommod_rest.services.interpro import interpro
from hommod_rest.services.blast import blaster
from hommod_rest.services.align import aligner

import domainalign

_log = logging.getLogger(__name__)


pFullyGapped = re.compile(r"^[\-\.]+$")
def is_fully_gapped(aligned_seq):
    return pFullyGapped.match(aligned_seq) is not None


# Record on how long certain modeling steps took,
# for performance monitoring:
def time_log(text):
    logfile = 'time.txt'
    open(logfile, 'a').write(text)


def select_highest_identity(seq, fromd):
    qid = 'query'

    bestID = 0
    best = None

    for homolog in fromd.keys():

        # Aligning them all in one run can be very time consuming !
        # So align in pairs

        aligned = aligner.clustal_align({qid: seq, homolog: fromd[homolog]})
        pid = getPercentageIdentity(aligned[qid], aligned[homolog])

        if pid > bestID:
            bestID = pid
            best = homolog

    return best


def write_alignment_fasta(chain_order, alignments_by_chain, template_ac, path):
    """
    Creates the input file for YASAYA's modeling run.
    Includes all the chain ids in 'chain_order' in the alignment.
    """
    f = open(path, 'w')

    f.write('>target\n')
    n = 0
    for chain in chain_order:
        if chain not in alignments_by_chain:
            continue
        if n > 0:
            f.write('|')
        f.write(alignments_by_chain[chain]['target'])
        n += 1
    f.write('\n')

    f.write('>%s\n' % template_ac)
    n = 0
    for chain in chain_order:
        if chain not in alignments_by_chain:
            continue
        if n > 0:
            f.write('|')
        f.write(alignments_by_chain[chain]['template'])
        n += 1
    f.write('\n')

    f.close()


# Modeling requires paths in the configuration settings.
# Therefore, required variables are wrapped in an object.
class Modeler(object):
    # must be set before building models:
    # * yasara dir: where yasara is installed
    # * execution root: where the modeler can run and create temporary files
    # * model root: where finished models should be placed
    def __init__(self, yasara_dir=None, execution_root_dir=None,
                 model_root_dir=None, template_blacklist=None):
        self._yasara_dir = yasara_dir
        self._execution_root_dir = execution_root_dir
        self._model_root_dir = model_root_dir
        self._template_blacklist = template_blacklist

    @property
    def yasara_dir(self):
        return self._yasara_dir

    @yasara_dir.setter
    def yasara_dir(self, yasara_dir):
        self._yasara_dir = yasara_dir

        if os.path.isdir(yasara_dir):
            # Need to set these path variables before importing the yasara module..
            sys.path.append(os.path.join(yasara_dir, 'pym'))
            sys.path.append(os.path.join(yasara_dir, 'plg'))

            # Don't use graphic mode, it fails without graphical interface!
            import yasaramodule
            self.yasara = yasaramodule
            self.yasara.info.mode = 'txt'

    @property
    def execution_root_dir(self):
        return self._execution_root_dir

    @execution_root_dir.setter
    def execution_root_dir(self, execution_root_dir):
        self._execution_root_dir = execution_root_dir

    @property
    def model_root_dir(self):
        return self._model_root_dir

    @model_root_dir.setter
    def model_root_dir(self, model_root_dir):
        self._model_root_dir = model_root_dir

    @property
    def template_blacklist(self):
        return self._template_blacklist

    @template_blacklist.setter
    def template_blacklist(self, template_blacklist):
        self._template_blacklist = template_blacklist

    def _check_init(self):
        if not self._yasara_dir:
            raise Exception("yasara_dir not set")

        if not self._execution_root_dir:
            raise Exception("execution_root_dir not set")

        if not self._model_root_dir:
            raise Exception("model_root_dir not set")

        if not self._template_blacklist:
            raise Exception("template_blacklist not set")

    # The blacklist contains pdb ids of templates that always
    # cause failures in yasara's modeling run.
    # pdbids are added to this blacklist file to exclude them
    # from future blast runs.
    def _add_template_to_blacklist(self, pdbid):

        self._check_init()

        _log.warn("adding template %s to blacklist %s" %
                  (pdbid.lower(), self.template_blacklist))

        with open(self.template_blacklist, 'a') as f:
            f.write('%s\n' % pdbid.lower())

    def _template_in_blacklist(self, pdbid):

        self._check_init()

        with open(self.template_blacklist, 'r') as f:
            pdbids = f.read().split()
            return pdbid.lower() in pdbids

    # Use this function to get the template sequences:
    def get_chain_order_and_sequences(self, tempobj):

        _log.info("indexing protein sequences present in yasara object")

        self._check_init()

        # chain_order is the order that yasara maintains.

        chain_order = self.yasara.ListMol('obj %i and protein' % tempobj, 'MOL')
        s = self.yasara.SequenceMol('obj %i and protein' % tempobj)

        uniquechain_order = []
        template_chain_sequences = {}
        for i in range(len(chain_order)):
            chain = chain_order[i]
            seq = ''

            for s in self.yasara.SequenceMol(
                    'obj %i and mol %s and protein' % (tempobj, chain)):
                seq += s

            # Some templates have two protein chains with the same ID:
            if chain not in template_chain_sequences:
                uniquechain_order.append(chain)

            template_chain_sequences[chain] = seq

        _log.info("found %d protein sequences present in yasara object"
                  % len(uniquechain_order))

        return [uniquechain_order, template_chain_sequences]

    # Get the sequence of one chain in YASARA.
    def getChainSeq(self, tempobj, chain):
        self._check_init()

        return self.yasara.SequenceMol('%s obj %i protein' %
                                       (chain, tempobj))[0]

    def _check_fix_chain_breaks(self, obj, chain):
        residues_atoms = {}
        res_order = []
        for s in self.yasara.ListAtom(
                        "obj %i and mol %s and aminoacid" % (obj, chain),
                        "RESNUM ATOMNAME"):

            # resnum is a string here, because it includes insertion codes!
            resnum, atomname = s.split()
            if resnum not in residues_atoms:
                residues_atoms[resnum] = []
                res_order.append(resnum)
            residues_atoms[resnum].append(atomname)

        prevresnum = None
        for resnum in res_order:
            incomplete = False
            for atomname in ['N', 'CA', 'C']:
                if atomname not in residues_atoms[resnum]:
                    _log.debug("residue {} of mol {} misses atom {}, it has only {}"
                               .format(resnum, chain, atomname, residues_atoms[resnum]))
                    incomplete = True
                    break

            # Delete residues with incomplete backbone:
            if incomplete:
                self.yasara.DelRes("obj %i and mol %s and res %s" % (obj, chain, resnum))
                continue

            if prevresnum is not None:
                # Connect all peptide bonds, even if N-C distance is huge.
                self.yasara.AddBond(
                    "obj %i mol %s res %s atom C" % (obj, chain, prevresnum),
                    "obj %i mol %s res %s atom N" % (obj, chain, resnum),
                    order=1.0
                )

            prevresnum = resnum


    # This function loads the template into yasara and also
    # takes some other actions:
    #    * oligomerization (can be switched off)
    #    * symmetry residues: fills in missing residues
    #    * deleting buffer molecules
    #    * replacing residues of type X by alanine or glycine
    #
    # returns the number of the yasara object and the
    # oligomerization factor. (1 if no oligomerization)
    def _set_template(self, tempac, oligomerize=True):
        _log.info("setting template %s to yasara" % tempac)
        self._check_init()

        self.yasara.Clear()

        tempobj = self.yasara.LoadPDB(tempac, download='yes')[0]
        self.yasara.DelObj('not %i' % tempobj)

        # Count the number of molecules in unoligomerized state:
        nMolsUnoligomerized = \
            len(self.yasara.ListMol('obj %i protein' % tempobj, 'MOL'))
        if nMolsUnoligomerized <= 0:
            raise Exception("No protein chains found in %s" % tempac)

        if oligomerize:
            try:
                # Oligomerisation might create extra objects,
                # if so join them together:
                oligolist = self.yasara.OligomerizeObj(tempobj)
                if len(oligolist) > 1:
                    for o in oligolist:
                        if o != tempobj:
                            self.yasara.JoinObj(o, tempobj)

                # oligomerisation = len(oligolist)
            except:
                # oligomerisation might throw an exception,
                # in which case we must restore the old situation
                self.yasara.Clear()
                tempobj = self.yasara.LoadPDB(tempac, download='yes')[0]

            # Count the number of molecules in oligomerized state:
            nMolsOligomerized = len(
                self.yasara.ListMol('obj %i protein' % tempobj, 'MOL'))

            # We want as many chains as possible in the template.
            # If oligomerisation decreased the number of molecules, undo it!
            if nMolsOligomerized < nMolsUnoligomerized:

                self.yasara.Clear()
                tempobj = self.yasara.LoadPDB(tempac, download='yes')[0]

        # If BuildSymRes throws an exception, we just continue
        try:
            self.yasara.BuildSymRes('obj %i' % tempobj)
        except Exception as e:
            _log.warn('Cannot execute BuildSymRes on {}: {}'.format(tempac,
                                                                    e.args[0]))
        # Make sure there's only one chain for each chain identifier:
        chain_order = self.yasara.ListMol('obj %i protein' % tempobj, 'MOL')
        for i in range(len(chain_order)):
            if i > 0 and chain_order[i - 1] == chain_order[i]:
                _log.debug("chain %s occurs twice in %s, joining.." %
                           (chain_order[i], tempac))

                atomnums = self.yasara.ListAtom('mol %s and obj %i and Protein' %
                                                (chain_order[i], tempobj), "ATOMNUM")
                # Call JoinMol on the last molecule in YASARA's order:
                self.yasara.JoinMol("atom %s" % atomnums[-1])

            self._check_fix_chain_breaks(tempobj, chain_order[i])

        # YASARA also cleans when starting a modeling run, but
        # we need to make sure that we have the molecule in its
        # final state, before we start reading from it.
        self.yasara.LogAs('clean.log')
        try:
            self.yasara.CleanObj(tempobj)
        except RuntimeError:
            # When yasara connection suddenly breaks without a reason,
            # look for errorexit.txt
            if os.path.isfile('errorexit.txt'):
                with open('errorexit.txt', 'r') as f:
                    raise Exception(f.read())

        # Some error messages can only be found in the log:
        with open('clean.log', 'r') as f:
            for line in f:
                if 'ERROR' in line:
                    raise Exception(line)

        # Make sure there are no chains with sequence XXXXXXXXXXXXXXXXXX,
        # Otherwise, yasara would remove the entire chain.
        self.yasara.SwapRes(
            'Protein and UNK and atom CA and atom CB', 'ALA')
        self.yasara.SwapRes(
            'Protein and UNK and atom CA and not atom CB', 'GLY')

        # Delete buffer molecules:
        self.yasara.DelRes('HOH H2O DOD D2O TIP WAT SOL ACT ACM ACY ' +
                           'EDO EOH FLC FMT GOL I3M IPA MLI MOH PEO ' +
                           'PO4 SO3 SO4 _TE UNX ACE')

        # Final checks for unfixable errors:
        chain_order = self.yasara.ListMol('obj %i protein' % tempobj, 'MOL')
        for chain in chain_order:
            nocc = 0
            for c in chain_order:
                if c == chain:
                    nocc += 1

            if nocc > 1:
                self._add_template_to_blacklist(tempac)
                seqs = self.yasara.SequenceMol('obj %i and mol %s' % (tempobj, chain))
                _log.error("chain {} in {} occurs more than once after cleaning: {}"
                           .format(chain, tempac, seqs))
                raise Exception("chain {} in {} occurs more than once after cleaning: {}"
                                .format(chain, tempac, seqs))

        _log.info("initialized yasara template with %d chains" % len(chain_order))

        return [tempobj, nMolsOligomerized / nMolsUnoligomerized]

    # This function must create the same file types as _build_for_domain,
    # but instead of building a model, use the pdb structure.
    def _collect_template(self, model_dir, main_target_id, uniprot_species_name,
                          main_template_id, main_domain_alignment, main_domain_range):
        self._check_init()
        _log.info("collecting template %s for %s, rather than building a model"
                  % (str(main_template_id), main_target_id))

        pdbac = main_template_id.pdbac.lower()
        chainID = main_template_id.chainID

        model_path = os.path.join(model_dir, 'target.pdb')
        alignmentFastaPath = os.path.join(model_dir, 'align.fasta')
        selectedTargetsPath = \
            os.path.join(model_dir, 'selected-targets.txt')
        try:
            # Download and decompress pdb file:
            open(model_path, 'w').write(get_pdb_contents(pdbac))

            # Document the match:
            open(selectedTargetsPath, 'w') \
                .write('template: %s_%s\n' % (pdbac, chainID))

            open(selectedTargetsPath, 'a') \
                .write('\tmatched main target %s with chain %s\n' %
                       (main_target_id, chainID))

            # Write down the alignment:
            f = open(alignmentFastaPath, 'w')
            f.write('>target\n')
            f.write(main_domain_alignment['target'] + '\n')
            f.write('>%s_%s\n' % (pdbac, chainID))
            f.write(main_domain_alignment['template'] + '\n')
            f.close()
        except:
            # Clean up files of failed operation:
            for path in [model_path, alignmentFastaPath, selectedTargetsPath]:
                if os.path.isfile(path):
                    os.remove(path)
            raise

    # This function builds a model for one alignment.
    # It must be given a directory to work in,
    # it also needs the main target's id, sequence and domain range object,
    # a uniprot species id, a template ID (pdbid and chain)
    def _build_for_domain(self, model_dir, main_target_id, uniprot_species_name,
                          main_template_id, main_target_sequence, main_domain_range):

        if main_domain_range.start >= main_domain_range.end:
            raise Exception("invalid main domain range: %d - %d" %
                    (main_domain_range.start, main_domain_range.end))
        if main_domain_range.start < 0 or \
                main_domain_range.end > len(main_target_sequence):
            raise Exception(
                    "main domain range: %d - %d exceeds sequence length(%d)" %
                    (main_domain_range.start, main_domain_range.end,
                     len(main_target_sequence)))

        self._check_init()
        _log.info("building model %s range %s on template %s" %
                  (main_template_id, str(main_domain_range), str(main_template_id)))

        model_path = os.path.join(model_dir, 'target.pdb')
        alignmentFastaPath = os.path.join(model_dir, 'align.fasta')
        selectedTargetsPath = \
            os.path.join(model_dir, 'selected-targets.txt')

        open(selectedTargetsPath, 'w') \
            .write('template: %s_%s\n' % (main_template_id.pdbac,
                                          main_template_id.chainID))

        # Load template and perform all necessary modifications:
        tempobj, oligomerisation = \
            self._set_template(main_template_id.pdbac)

        # One alignment per chain id:
        alignments = {}

        # Get the template's sequences:
        chain_order, template_chain_sequences = \
            self.get_chain_order_and_sequences(tempobj)

        # Remove funny insertions/deletions from main target:
        main_target_sequence = \
            adjust_target(main_target_sequence,
                          template_chain_sequences[main_template_id.chainID],
                          uniprot_species_name)

        _log.debug("template sequences after _set_template:\n " +
                   str(template_chain_sequences))

        # A template can have more than one chain, but the main target is
        # not necessarily homologous to all of these chains. So we must make
        # sure that we find other targets for those other template chains.
        # However, we only need to do this for chains that actually interact
        # with the main target's template chain...

        time_start = time()

        # Part of the main target sequence that's in the domain range:
        mainDomainSeq = \
            main_target_sequence[main_domain_range.start: main_domain_range.end]

        # Add the first alignments,
        # involving the main target sequence:
        # (the main target might be selected for more than one template sequence)
        main_target_chain_ids = pick_template_chains_for(template_chain_sequences, main_target_sequence, main_template_id)
        for chainID in main_target_chain_ids:

            # Of coarse, we've already made an alignment previously
            # while determining the covered domain. However, we want
            # to make a new final alignment here, using a fine-tuned
            # program.

            # It's based on secondary structure, so we get that
            # information first:
            tempCAs, templateChainSeq, templateChainSecStr = \
                getChainCAsSeqSecStr(self.yasara, tempobj, chainID)

            _log.debug("aligning input target sequence \n{}"
                       .format(mainDomainSeq) +
                       "\nto template {} {}"
                       .format(main_template_id.pdbac, chainID) +
                       " sequence\n{}\n{}"
                       .format(templateChainSeq, templateChainSecStr))

            alignments[chainID] = aligner \
                .kmad_align(templateChainSeq,
                            templateChainSecStr, mainDomainSeq)
            if is_fully_gapped(alignments[chainID]['target']):
                raise Exception(
                         "fully gapped main target alignment returned by kmad")

            open(selectedTargetsPath, 'a') \
                .write('\tmodeling main target %s on chain %s\n' %
                       (main_target_id, chainID))

        # Try to find and align target sequences for
        # interacting chains in the template, while keeping
        # in mind which residues interact and must thus be
        # covered by the alignment.
        # We expand the set of involved template
        # chains with every iteration, until all template
        # chains have been added.
        yasara_chains = {}
        for c in self.yasara.ListMol("obj %i and protein" % tempobj, "MOL"):
            yasara_chains[c] = YasaraChain(self.yasara, tempobj, c)
        while len(alignments) < \
                len(self.yasara.ListMol('obj %i and protein' % tempobj)):

            # First, make python remember to which
            # chains the candidate chains interact:
            candidateChainInteractsWith = {}
            for c in alignments:
                for chainID in \
                        listInteractingChains(yasara_chains[c]):

                    # Skip those, that we've already aligned
                    # to prevent infinite loops:
                    if chainID in alignments:
                        continue

                    if chainID not in candidateChainInteractsWith:
                        candidateChainInteractsWith[chainID] = []

                    candidateChainInteractsWith[chainID].append(c)

            if len(candidateChainInteractsWith) <= 0:
                _log.debug("No more interacting chains to add, quitting loop")
                break

            # iterate over chains that might interact with the chains
            # that are already in the set:
            for chainID in candidateChainInteractsWith:

                # Gather the alignments of
                # the interaction partner chains.
                # These are chains that we've already added to the set,
                # so their alignments are already there.
                interacting_chain_alignments = {}
                for interactingChainID in \
                        candidateChainInteractsWith[chainID]:

                    interacting_chain_alignments[interactingChainID] = \
                        alignments[interactingChainID]

                tempCAs, templateChainSeq, templateChainSecStr = \
                    getChainCAsSeqSecStr(self.yasara, tempobj, chainID)

                # Find targets from the target species, homologous to this template chain:
                potentialtarget_seqs = find_targets(templateChainSeq,
                                                    uniprot_species_name)

                _log.debug("found {} potential target sequences for {} {}"
                           .format(len(potentialtarget_seqs),
                                   main_template_id.pdbac, chainID))

                yasaraChain = yasara_chains[chainID]

                # Pick the target with the hightest sequence identity,
                # and domain coverage.
                targetsInterproRanges = {}
                bestPID = 0.0
                selectedTarget = None
                for targetID in potentialtarget_seqs.keys():

                    _log.debug(
                        "aligning potential target sequence\n{}\n"
                        .format(potentialtarget_seqs[targetID]) +
                        " to interacting template {} {}"
                        .format(main_template_id.pdbac, chainID) +
                        " sequence\n{}\n{}"
                        .format(templateChainSeq, templateChainSecStr))

                    # If we have targets with extremely high coverage,
                    # then we don't need interpro
                    alignment = aligner.kmad_align(
                        templateChainSeq,
                        templateChainSecStr,
                        potentialtarget_seqs[targetID])

                    pcov, pid = getCoverageIdentity(
                        alignment['template'], alignment['target'])

                    if pcov > 90.0:
                        # Coverage is large enough, take the whole sequence:

                        if chainID not in alignments or pid > bestPID:
                            alignments[chainID] = alignment
                            bestPID = pid
                            selectedTarget = targetID

                    else:
                        # no luck, need to bother interpro:
                        targetsInterproRanges[targetID] = interpro \
                            .get_domain_locations(main_target_sequence)

                if chainID not in alignments:  # Not been added yet at this point

                    _log.debug("no 90 percent coverage target for {} {}"
                               .format(main_template_id.pdbac, chainID))

                    # Find out which targets, given their domain ranges,
                    # preserve the interaction with the other chains when
                    # aligned to the template chain:
                    picker = InteractionPicker(chainID,
                                               yasara_chains,
                                               interacting_chain_alignments)
                    alignmentTriplesPerTarget = \
                        domainalign.pickAlignments(
                            yasaraChain, potentialtarget_seqs,
                            targetsInterproRanges, picker)

                    _log.debug("got {} targets passed by domainalign"
                               .format(len(alignmentTriplesPerTarget)) +
                               " interacton picker for {} {}"
                               .format(main_template_id.pdbac, chainID))

                    # Iterate over targets that passed:
                    for targetID in alignmentTriplesPerTarget:

                        nrange = len(alignmentTriplesPerTarget[targetID])
                        _log.debug("got {} interpro ranges for {} on {} {}"
                                   .format(nrange, targetID,
                                           main_template_id.pdbac, chainID))

                        # Join the alignments for this one target into one,
                        # maximizing the coverage of the template chain:
                        domain_alignment = domainalign \
                            .joinAlignmentsToBestTemplateCoverage(
                                alignmentTriplesPerTarget[targetID])

                        # Remove gaps from aligned target sequence:
                        domain_target_seq = \
                            domain_alignment['target'].replace('-', '')

                        _log.debug(
                            "aligning {} domain target sequence\n{}"
                            .format(targetID, domain_target_seq) +
                            "\nto interacting template {} {}"
                            .format(main_template_id.pdbac, chainID) +
                            " sequence\n{}\n{}"
                            .format(templateChainSeq, templateChainSecStr))

                        # Realign once more using the complete target
                        # range that we found by joining the alignments:
                        alignment = aligner.kmad_align(
                            templateChainSeq, templateChainSecStr,
                            domain_target_seq)

                        # Determine the best target for this chain id:
                        # (highest identity)
                        nalign, pid = getNalignIdentity(
                            alignment['target'], alignment['template'])
                        if chainID not in alignments or pid > bestPID:
                            alignments[chainID] = alignment
                            bestPID = pid
                            selectedTarget = targetID

                _log.debug("selected target for interacting {} {} is {}"
                           .format(main_template_id.pdbac, chainID,
                                   selectedTarget))

                # Occasionally, we don't find any target for a template chain.
                # It's solved here:
                if chainID not in alignments:

                    _log.debug("putting poly-A on chain {} of {}, sequence=\n\"{}\""
                               .format(chainID, main_template_id.pdbac, templateChainSeq))

                    # Place alanines in the target sequence, with same
                    # length as templare chain:
                    alignments[chainID] = {
                        'template': templateChainSeq,
                        'target': 'A' * len(templateChainSeq)}
                    selectedTarget = 'poly-A'
                    self.yasara.OccupMol(chainID, 0.0)

                # We document in a file which targets we used for this model:
                open(selectedTargetsPath, 'a') \
                    .write('\tmodeling target %s on chain %s\n'
                           % (selectedTarget, chainID))

            # < end of interaction finding for-loop

        # Delete chains that weren't aligned, assuming there's no
        # interaction with the main target's homologs:
        for chainID in chain_order:
            if chainID in alignments and \
                    is_fully_gapped(alignments[chainID]['target']):
                alignments.pop(chainID)

            if chainID not in alignments:

                _log.debug("deleting not-interacting chain {} of {}"
                           .format(chainID, main_template_id.pdbac))

                # Only delete protein, keep ligands and nucleic acids, for
                # which no alignments are created.
                self.yasara.DelMol('%s and protein' % chainID)
                open(selectedTargetsPath, 'a') \
                    .write('\tdeleting not-interacting ' +
                           'chain %s from template\n' % chainID)

        while True:
            # The set of chains might have changed:
            chain_order, template_chain_sequences = \
                self.get_chain_order_and_sequences(tempobj)

            _log.debug("Writing down alignments for chain order: {}".format(chain_order))

            # Make the alignment file for yasara:
            write_alignment_fasta(
                chain_order, alignments, main_template_id.pdbac,
                alignmentFastaPath)

            # Start the modeling run:
            self.model_with_alignment(alignmentFastaPath, tempobj)

            if os.path.isfile("target.yob"):

                # Save the model in PDB format:
                self.yasara.SavePDB(tempobj, model_path)
                _log.info("sucessfully created " + model_path)

                return
            else:  # target.yob is missing

                new_chain_order, new_template_chain_sequences = \
                    self.get_chain_order_and_sequences(tempobj)

                if len(new_chain_order) <= 0:
                    raise RuntimeError("yasara experiment removed all chains")

                elif new_chain_order != chain_order or new_template_chain_sequences != template_chain_sequences:
                    _log.debug("the chains have been modified during the experiment, retrying..")

                    # Redo the alignments with the new chains:
                    for chainID in new_chain_order:
                        tempCAs, templateChainSeq, templateChainSecStr = \
                            getChainCAsSeqSecStr(self.yasara, tempobj, chainID)

                        alignments[chainID] = aligner.kmad_align(
                            templateChainSeq, templateChainSecStr,
                            alignments[chainID]['target'])

                    continue

                elif os.path.isfile("errorexit.txt"):
                    with open("errorexit.txt", 'r') as f:
                        raise RuntimeError(f.read())
                else:
                    raise RuntimeError("yasara modeling run did not complete for %s %s (%d - %d)\n%s\n\n"
                                       % (uniprot_species_name, main_template_id,
                                          main_domain_range.start, main_domain_range.end,
                                          main_target_sequence) +
                                       "\'target.yob\' is missing\n" +
                                       "Please check yasara's output for further details")

    def modelProc(self, main_target_sequence, uniprot_species_name, requireRes=None,
                  overwrite=False, chosenTemplateID=None):
        """
        This is the main function for building a set of models.
        Input is a target sequence and a uniprot species id.
        The function searches its own templates.

        Optional:
            requireRes: a residue number in target sequence (1,2,3,..) that
            must be covered by the model.
            overwrite: true if old models must be rebuilt, false if they must
            be skipped
        """
        self._check_init()

        _log.info("building models for %s in %s" % (main_target_sequence, uniprot_species_name))

        if not os.path.isdir(self.model_root_dir):
            os.mkdir(self.model_root_dir)

        if not os.path.isdir(self.execution_root_dir):
            os.mkdir(self.execution_root_dir)

        # A model can have multiple target sequences, since the template
        # can have multiple chains.
        # The target sequence, that the model was built for, will be
        # called the 'main' target sequence.

        # Unique ID for this sequence:
        main_target_id = idForSeq(main_target_sequence)

        # Determine all domain ranges within our main target sequence.
        # (source: interpro)
        ranges = interpro.get_domain_locations(main_target_sequence)

        # yasara sticks to the directory where it was started,
        # so make a unique directory for yasara to run in and
        # let it store all its output files there:
        run_dir = os.path.join(self.execution_root_dir,
                               'run-yasara-%s-%s-%i' % (uniprot_species_name,
                                                        main_target_id,
                                                        os.getpid()))
        if os.path.isdir(run_dir):
            _log.debug("removing old verion of %s" % run_dir)
            shutil.rmtree(run_dir)
        os.mkdir(run_dir)
        os.chdir(run_dir)

        # Restart yasara, in case if it was already running.
        # We don't want it running in a different directory,
        # could cause errors.
        if self.yasara.pid:
            _log.debug("restarting yasara for %s %s" % (uniprot_species_name,
                                                        main_target_sequence))
            self.yasara.Exit()
            self.yasara.start()

        # Blast and filter alignments that satisfy the given set of
        # domain ranges. No halfway cut domains are alowed!
        if chosenTemplateID:
            tempobj, oligomerisation = \
                self._set_template(chosenTemplateID.pdbac)
            yasaraChain = YasaraChain(self.yasara, tempobj,
                                      chosenTemplateID.chainID)
            maintarget_alignments = \
                domainalign.getAlignments(ranges, main_target_sequence, yasaraChain)
        else:
            maintarget_alignments = \
                domainalign.getAlignments(ranges, main_target_sequence)

        if len(maintarget_alignments) <= 0:
            _log.warn('no alignments found for sequence:\n' + main_target_sequence)

        _log.info('got %d alignments for sequence:\n%s'
                  % (len(maintarget_alignments), main_target_sequence))

        model_paths = []
        failed_models = {}

        realign = False

        # Iterate over all alignments that we've got. Any alignment
        # is a potential model.
        for main_domain_range, main_template_id, main_domain_alignment in \
                maintarget_alignments:

            # Make sure that we're still in the run directory after the previous iteration:
            os.chdir(run_dir)

            # skip ranges that don't cover the requireRes (if given)
            if requireRes and \
                    ((requireRes - 1) < main_domain_range.start or
                     (requireRes - 1) >= main_domain_range.end):

                _log.debug("skipping range %d - %d on %s, bause it does not cover residue %d\n%s\n%s"
                           % (main_domain_range.start, main_domain_range.end, main_template_id, requireRes,
                              main_domain_alignment['target'], main_domain_alignment['template']))
                continue

            if chosenTemplateID:
                model_name = '%s_%s_%i-%i_%s' % \
                    (main_target_id, uniprot_species_name,
                     main_domain_range.start + 1, main_domain_range.end,
                     str(chosenTemplateID))
            else:
                model_name = '%s_%s_%i-%i' % \
                    (main_target_id, uniprot_species_name,
                     main_domain_range.start + 1, main_domain_range.end)

            _log.debug("locking %s, running from %s" % (model_name, run_dir))

            model_dir = os.path.join(self.model_root_dir, model_name)
            model_archive = model_dir + '.tgz'

            # A unique lockfile name, based on the model name:
            lockfile_path = model_dir + '_lock'

            # we dont want two threads building the same model, so lock:
            lock = FileLock(lockfile_path)
            with lock:
                _log.debug("locked {}".format(lockfile_path))
                if os.path.isfile(model_archive) and not overwrite:

                    # Found the archive for this model before starting,
                    # see if it's populated:

                    members = tarfile.open(model_archive, mode='r:gz').getnames()
                    pdbpath = os.path.join(model_name, 'target.pdb')
                    if pdbpath in members:

                        model_paths.append(model_archive)

                        _log.info('%s exists, skipping..' % model_archive)
                        continue

                # the model directory is going to be turned into an archive,
                # All created files will be sent there eventually,
                # but it's not yasara's working directory.
                if not os.path.isdir(model_dir):
                    os.mkdir(model_dir)

                model_path = os.path.join(model_dir, 'target.pdb')

                if os.path.isfile(model_path) and not overwrite:
                    # Found the output directory for this model already
                    # before we started. (just not yet archived)

                    # archive and clean up:
                    os.chdir(self.model_root_dir)
                    tf = tarfile.open(model_name + '.tgz', 'w:gz')
                    tf.add(model_name)  # refers to model_dir
                    tf.close()
                    shutil.rmtree(model_dir)

                    model_paths.append(model_archive)
                    _log.info('%s already exists, skipping..' % model_archive)
                    continue

                # We're now sure that the model doesn't exist yet.
                # At this point we start the actual model building.
                try:
                    # If the complete template sequence fully matches the target, then
                    # there's no need to build the model.

                    nalign, pid = getNalignIdentity(main_domain_alignment['target'],
                                                    main_domain_alignment['template'])
                    template_seq = getTemplateSequence(main_template_id.pdbac,
                                                       main_template_id.chainID)
                    if pid >= 100.0 and template_seq == main_domain_alignment['template']:

                        _log.debug("now calling _collect_template for %s %s (%d - %d)\n>target\n%s\n>template\n%s"
                                   % (uniprot_species_name, main_template_id,
                                      main_domain_range.start, main_domain_range.end,
                                      main_domain_alignment['target'], main_domain_alignment['template']))

                        self._collect_template(model_dir, main_target_id,
                                               uniprot_species_name, main_template_id,
                                               main_domain_alignment, main_domain_range)

                    else:
                        _log.debug("now calling _build_for_domain for %s %s (%d - %d)\n%s"
                                   % (uniprot_species_name, main_template_id,
                                      main_domain_range.start, main_domain_range.end,
                                      main_target_sequence))

                        self._build_for_domain(model_dir, main_target_id,
                                               uniprot_species_name, main_template_id,
                                               main_target_sequence, main_domain_range)

                except Exception as ex:  # An error ocurred during modeling.

                    # Don't exit, just print error to log and
                    # move on to the next alignment.

                    exc_type, exc_value, exc_traceback = sys.exc_info()
                    stacktrace = ''.join(traceback.format_exception(
                        exc_type, exc_value, exc_traceback))

                    _log.error('an exception occured for {}:\n{}'
                               .format(model_name, stacktrace))

                    if self._template_in_blacklist(main_template_id.pdbac):

                        # The template has just been blacklisted, need to find a new one
                        realign = True
                    else:
                        failed_models[model_name] = stacktrace

                    # Also include error and yasara scene in the model dir,
                    # for debugging and yasara bug reports.
                    open(os.path.join(model_dir, 'errorexit.txt'), 'w') \
                        .write(stacktrace)
                    try:
                        self.yasara.SaveSce(os.path.join(model_dir, 'errorexit.sce'))
                    except:
                        # connection with yasara isn't alive anymore !
                        raise ex

                # At this point the model building has finished.
                # Move files from temporary run directory to their final destination..

                # Move all the files that yasara created:
                for f in os.listdir(run_dir):
                    os.rename(os.path.join(run_dir, f),
                              os.path.join(model_dir, f))

                # Create archive:
                parent = os.path.dirname(model_dir)
                os.chdir(parent)
                tf = tarfile.open(model_name + '.tgz', 'w:gz')
                tf.add(model_name)  # refers to model_dir
                tf.close()
                shutil.rmtree(model_dir)
                model_paths.append(model_archive)
                os.chdir(run_dir)

            # end of this iteration, move over to next range ...
            _log.debug("ending lock on {}".format(lockfile_path))

        # Clean up all runtime files:
        _log.debug("modeling done, cleaning up %s" % run_dir)
        os.chdir(self.execution_root_dir)
        shutil.rmtree(run_dir)

        _log.info("finished building all models for %s in %s" %
                  (main_target_sequence, uniprot_species_name))

        if realign:
            # Repeat once more with different templates:
            return self.modelProc(main_target_sequence, uniprot_species_name,
                                  requireRes, overwrite, chosenTemplateID)

        elif len(failed_models) > 0:
            s = ''
            for name in failed_models:
                s += "%s:\n%s\n" % (name, failed_models[name])

            raise Exception("the following models have failed:\n" + s)
        else:
            return model_paths

    def check_model(self, model_obj, template_ac, expected_chains, target_alignment,
                    template_alignment):
        self._check_init()
        """
        This function is called after yasara has finished builing the model.
        It must perform the final checks to be sure that the model is not too
        bad. If everything is OK, it must return nothing. If something is
        wrong, it returns a string telling what's wrong.
        """

        _log.info("checking finished model")

        # Check that all chains are there:
        chain_order, template_chain_sequences = self.get_chain_order_and_sequences(model_obj)
        if len(chain_order) < len(expected_chains):
            return "missing %i chains" % (len(expected_chains) - len(chain_order))

        # Verify that the alignment is good enough:
        # include one extra iteration to finalize the last chain
        nalign = 0
        nid = 0
        nseq = 1
        # include one extra iteration to finalize the last chain
        for i in range(len(target_alignment) + 1):

            if i >= len(target_alignment) or target_alignment[i] == '|':
                # end of chain

                pid = (100.0 * nid) / nalign
                if pid < min(95, minIdentity(nalign)):
                    return ("sequence %i identity(%.1f)" +
                            " too low for %i residues") % (nseq, pid, nalign)

                nseq += 1
                nalign = 0
                nid = 0
                continue

            else:
                if target_alignment[i].isalpha() and \
                        template_alignment[i].isalpha():
                    nalign += 1
                    if target_alignment[i].upper() == \
                            template_alignment[i].upper():
                        nid += 1

        # Verify that packing quality didn't decrease to much,
        # relative to template:
        pdbname = 'model%i.pdb' % os.getpid()
        self.yasara.SavePDB(model_obj, pdbname)

        whatifClient = client.Client('http://wiws.cmbi.ru.nl/wsdl',
                                     timeout=60 * 60)

        # Get packing qualities from whatif and determine how much
        # worse the model is than the template:
        modelWhatifID = whatifClient.service.UploadPDB(
            open(pdbname, 'r').read())

        os.remove(pdbname)

        qr = whatifClient.service.PackingQualityMolecule(modelWhatifID)
        if hasattr(qr[0], 'value'):
            modelPackingQuality = qr[0].value
            open('packing.log', 'a') \
                .write('model packing quality is %.3f\n' % modelPackingQuality)
        else:

            _log.warn("unexpected result from whatif:\n" + str(qr))
            return 'whatif output for model: ' + str(qr)

        qr = whatifClient.service.PackingQualityMolecule(template_ac)
        if hasattr(qr[0], 'value'):
            templatePackingQuality = qr[0].value
            open('packing.log', 'a') \
                .write('template packing quality is %.3f\n' %
                       templatePackingQuality)
        else:
            _log.warn("unexpected result from whatif:\n" + str(qr))
            return 'whatif output for template: ' + str(qr)

        if modelPackingQuality < (templatePackingQuality - 1.0):

            _log.warn("Model packing quality %.3f too low" +
                      ", compared to %.3f from template" %
                      (modelPackingQuality, templatePackingQuality))

            return ("Model packing quality %.3f too low" +
                    ", compared to %.3f from template" %
                    (modelPackingQuality, templatePackingQuality))

    # This function starts the model building process by calling
    # yasara. It will wait for the modeling run to finish.
    def model_with_alignment(self, alignmentFastaPath, tempobj):
        _log.info("building a model in yasara from alignment %s" % alignmentFastaPath)
        _log.debug("running yasara homology modeling experiment from " + os.getcwd())

        if not os.path.isdir(os.getcwd()):
            raise Exception("current work dir has been deleted")

        self._check_init()

        self.yasara.Processors(1)

        self.yasara.LogAs('model.log')
        self.yasara.ExperimentHomologyModeling(
            templateobj=tempobj,
            alignfile=alignmentFastaPath,
            templates="1, sameseq = 1",
            alignments=1,
            termextension=0,
            oligostate=32,
            looplenmax=10,
            animation='fast',
            speed='fast',
            loopsamples=20,
            resultfile='target'
        )
        self.yasara.Experiment("On")
        self.yasara.Wait("Expend")

        with open('model.log', 'r') as log:
            for line in log:
                if 'ERROR' in line:
                    raise Exception(line)

        _log.info("yasara modeling run ended for alignment %s"
                  % alignmentFastaPath)


modeler = Modeler()


def get_isoforms(seq, uniprot_species_name):
    """
    This function finds uniprot isoforms for a given sequence.
    (a sequence that looks almost like it, except for some deletions)
    It makes sure the isoform comes from the given species.
    """

    hits = blaster.blast_species(seq, uniprot_species_name)

    selected = {}
    for hitID in hits:
        for ali in hits[hitID]:
            ac = hitID.split('|')[1]
            if ali.getIdentity() >= 99.0:
                if ac not in selected or \
                        selected[ac].getNumberResiduesAligned() < \
                        ali.getNumberResiduesAligned():
                    selected[ac] = ali

    _log.info("found %s isofroms %s for\n%s" % (uniprot_species_name,
                                                str(selected), seq))

    return selected


def adjust_target(target_seq, templateseq, uniprot_species_name):
    """
    This function fills in deletions with pieces from isoforms in uniprot.
    It's necessary sometimes, because deletions might destroy the
    hydrofobic core of a model.
    """

    seqs = {}
    tarkey = 'target'
    seqs[tarkey] = target_seq
    tempkey = 'template'
    seqs[tempkey] = templateseq

    isoformkeys = []
    isoformHits = get_isoforms(target_seq, uniprot_species_name)
    for ac in isoformHits.keys():
        hit = isoformHits[ac]
        seqs[ac] = hit.subjectalignment.replace('-', '').upper()
        if len(seqs[ac]) < 10000:  # avoid titin
            isoformkeys.append(ac)

    aligned = aligner.clustal_align({
        tempkey: seqs[tempkey], tarkey: seqs[tarkey]})
    deletion_ranges = identifyDeletedRegions(aligned[tarkey])

    if len(deletion_ranges) <= 0:
        return target_seq

    adjusted = aligned[tarkey]

    # insertions must be added in reverse order, so that indices don't get
    # messed up
    deletion_ranges.reverse()
    for i, f in deletion_ranges:

        # i:f marks the start and end index in the alignment
        # between template and target

        if len(aligned[tempkey][i:f].replace('-', '')) == 0:

            # Template region is empty too, not a real deletion,
            # no alignment possible
            continue

        # deletionLength = f - i

        tarn = aligned[tarkey][:i].replace('-', '')
        tarc = aligned[tarkey][f:].replace('-', '')

        naligns = {}  # counts residues aligned to the template per isoform
        pids = {}  # % identity to template for each isoform
        best = None  # Best replacement key

        replacement = aligned[tarkey][i:f]  # Initialize as a gap

        for key in isoformkeys:

            # Compare the target sequence to the isoform
            # to see if pieces can be added:

            alignedisoform = \
                aligner.clustal_align({tarkey: seqs[tarkey], key: seqs[key]})

            # identify the location of the deletion in the isoform
            ai = 0
            while alignedisoform[tarkey][:ai].replace('-', '') != tarn:
                ai += 1

            af = len(alignedisoform[tarkey]) - 1
            while alignedisoform[tarkey][af:].replace('-', '') != tarc:
                af -= 1

            # ai:af marks the start and end index in the alignment
            # between template and isoform

            naligns[key] = 0
            nid = 0
            if len(alignedisoform[key][ai:af].replace('-', '')) > 0:

                # Align to template fragment:
                ali = aligner.clustal_align({
                    key: alignedisoform[key][ai:af].replace('-', ''),
                    tempkey: aligned[tempkey][i:f].replace('-', '')
                })

                for a in range(len(ali[key])):
                    if ali[key].isalpha() and ali[tempkey].isalpha():
                        naligns[key] += 1
                        if ali[key][a].upper() == ali[tempkey][a].upper():
                            nid += 1

            # Calculate sequence identity to template:
            if naligns[key] > 0:
                pids[key] = (100.0 * nid) / naligns[key]
            else:
                pids[key] = 0.0

            # the more residues filled in, the better
            # second most important is % identity
            if not best or naligns[key] > naligns[best] or \
                    (naligns[key] == naligns[best] and pids[key] > pids[best]):
                best = key
                replacement = alignedisoform[best][ai:af]

        # Fill the deletion with the chosen replacement piece:
        adjusted = adjusted[:i] + replacement + adjusted[f:]

    _log.info("adjusted\n%s\nto\n%s" % (aligned[tarkey], adjusted))

    return adjusted.replace('-', '')


# For a given sequence, finds an ortholog in the requested species:
# Criteria to be an ortholog: >70% identity and >90% coverage
def find_orthologs_to_sequence(seq, uniprot_species_name):
    orthologs = {}

    hits = blaster.blast_species(seq, uniprot_species_name)
    for hitID in hits.keys():
        for alignment in hits[hitID]:

            pcov, pid = getCoverageIdentity(
                alignment.queryalignment, alignment.subjectalignment)

            if pid > 70.0:

                ac = hitID.split('|')[1]
                oseq = getUniprotSeq(ac)

                if pcov > 90.0:
                    orthologs[ac] = oseq

    _log.info("found {} {} orthologs to\n{}"
              .format(len(orthologs), uniprot_species_name, seq))

    return orthologs


# The following function is responsible for making the first selection of
# target sequences tht could be modeled on a given template chain sequence.
# The selection will be refined later, based on domain-related criteria.
def find_targets(template_chain_sequence, uniprot_species_name):
    return find_orthologs_to_sequence(template_chain_sequence,
                                      uniprot_species_name)


# This function is responsible for grouping sequences that are more than
# 99% identical. The input is a dictionary of sequences and the
# output is a nested list of that dictionary's keys.
def group_identicals(d):
    if len(d) <= 1:
        return d.keys()

    alignments = {}

    grouped = []

    ids = d.keys()
    while len(ids) > 0:
        ID = ids[0]
        ids.remove(ID)
        grouped.append([ID])

        for otherID in ids[:]:

            pair = [ID, otherID]
            pair.sort()
            pair = pair[0] + '_' + pair[1]

            # Aligning them all in one run can be very time-consuming

            if pair not in alignments:
                alignments[pair] = aligner.clustal_align({
                    ID: d[ID],
                    otherID: d[otherID]
                })
            alignedD = alignments[pair]

            if getPercentageIdentity(alignedD[ID], alignedD[otherID]) >= 99.0:
                grouped[-1].append(otherID)
                ids.remove(otherID)

    _log.info("grouped %d sequences to %d groups" % (len(d), len(grouped)))

    return grouped


def pick_template_chains_for(chain_sequences, target_seq, chosen_template_id):
    """
    This function is used for choosing which chains of the template will be
    used to model the main target sequence on.
    It's given a main target sequence and all the chains of the template.
    It will return a list with one chain id or more. (in case of duplicates)
    """

    identicalChains = group_identicals(chain_sequences)
    hit_chain = chosen_template_id.chainID

    for group in identicalChains:
        if hit_chain in group:
            _log.info("selected chain %s in group %s for\n%s"
                      % (hit_chain, len(group), target_seq))
            return group

    _log.error('no hit found for %s in template' % target_seq)
    raise Exception('no hit found for %s in template' % target_seq)
