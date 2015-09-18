import logging
import os
import re
import shutil
import sys
import tarfile
import traceback

from time import time
from suds import client

from modelutils import (getChainCAsSeqSecStr, getNalignIdentity,
                        getCoverageIdentity, YasaraChain,
                        identifyDeletedRegions, getPercentageIdentity,
                        filterMinLength, minIdentity,
                        getUniprotSeq, idForSeq)

from interaction import listInteractingChains, InteractionPicker

from filelock import FileLock
from hommod_rest.services.interpro import interpro
from hommod_rest.services.blast import blaster
from hommod_rest.services.align import aligner

import domainalign

_log = logging.getLogger(__name__)

sh = logging.StreamHandler()
formatter = logging.Formatter(
    '%(asctime)s - %(levelname)s - %(message)s')
sh.setFormatter(formatter)
_log.addHandler(sh)
_log.setLevel(logging.DEBUG)


def time_log (text):

    logfile = 'time.txt'
    open (logfile, 'a').write (text)

        
def selectHighestSeqID(seq, fromd):
    qid = 'query'

    bestID = 0
    best = None

    for homolog in fromd.keys():

        # Aligning them all in one run can be very time consuming !
        # So align in pairs

        aligned = aligner.clustalAlign({qid: seq, homolog: fromd[homolog]})
        pid = getPercentageIdentity(aligned[qid], aligned[homolog])

        if pid > bestID:
            bestID = pid
            best = homolog

    return best


# Returns true if a sequence is only X's
def onlyX (seq):
    return (len(seq) == seq.count('X'))


# Many templates contain duplicate chains. This function reduces
# The set of sequences to a smaller set with only unique sequences.
def withoutHomomerDuplicates (chainOrder, chainSequences):

    count = {}
    seq2chains = {}

    # Order chain ids by sequence:
    for chain in chainOrder:
        seq = chainSequences [chain]
        if seq not in count:
            count [seq] = 0
            seq2chains [seq] = []
        count [seq] += 1
        seq2chains [seq].append (chain)

    # Check that every sequence occurs the same number of times:
    sameNum = True
    for n in count.values()[1:]:
        if n != count.values()[0]:
            sameNum = False
            break

    if sameNum: # If all sequences occur twice for example:

        # Return only the first of every set of duplicates:
        newChainSequences = {}
        newChainOrder = []
        for seq in count.keys ():
            newChainSequences [seq2chains [seq][0]] = seq

        for chain in chainOrder:
            if chain in newChainSequences:
                newChainOrder.append (chain)

        return [newChainOrder, newChainSequences]
    else:
        return [chainOrder, chainSequences]


def parseAli(filename, tempac):
    targetp = re.compile('\[target(.+)?\]$')

    seqs = {}
    ID = None
    for line in open(filename, 'r').readlines():
        if line.startswith('Sequence:'):

            line = line.strip()
            m = targetp.search(line.strip())
            if m:
                seqs['target'] = line[:m.start()].split(':')[1].strip()
            elif ID.upper().startswith(tempac.upper() + '-'):
                seqs[tempac] = line.split(':')[1].strip()[:len(seqs['target'])]
        else:
            ID = line.split(':')[0].strip()

    return seqs


def writeChainOrigins(description, filename):
    f = open(filename, 'w')
    for chain in description.keys():
        f.write('%s:%s\n' % (chain, description[chain]))
    f.close()


def writeDomainMap(filepath, targetAlign, templateAlign, ok):
    f = open(filepath, 'w')
    for i in range(0, len(targetAlign), 100):
        j = min(i + 100, len(targetAlign))
        f.write(targetAlign[i:j] + '\n')
        for n in ok[i:j]:
            if n < 0:
                f.write('-')
            elif n == 0:
                f.write('?')
            elif n > 0:
                f.write('+')
            else:
                f.write(' ')
        f.write('\n'+templateAlign[i:j]+'\n\n')
    f.close()


def writeAlignmentFasta(chainOrder, alignmentsByChain, templateAC, path):
    f = open(path, 'w')

    f.write('>target\n')
    n = 0
    for chain in chainOrder:
        if n > 0:
            f.write('|')
        f.write(alignmentsByChain[chain]['target'])
        n += 1
    f.write('\n')

    f.write('>%s\n' % templateAC)
    n = 0
    for chain in chainOrder:
        if n > 0:
            f.write('|')
        f.write(alignmentsByChain[chain]['template'])
        n += 1
    f.write('\n')

    f.close()


class Modeler(object):

    # must be set before building models:
    # * yasara dir: where yasara is installed
    # * execution root: where the modeler can run and create temporary files
    # * model root: where finished models should be placed

    def __init__(self, yasara_dir=None):
        self._yasara_dir = yasara_dir
        self.execution_root_dir = None
        self.model_root_dir = None
        self.template_blacklist = None

    @property
    def yasara_dir(self):
        return self._yasara_dir

    @yasara_dir.setter
    def yasara_dir(self, yasara_dir):
        self._yasara_dir = yasara_dir

        if not os.path.isdir(yasara_dir):
            raise ValueError("{} not found".format(yasara_dir))

        sys.path.append(os.path.join(yasara_dir, 'pym'))
        sys.path.append(os.path.join(yasara_dir, 'plg'))

        import yasaramodule
        self.yasara = yasaramodule
        self.yasara.info.mode = 'txt' # Don't use graphic, it fails without graphical interface

    def _check_init(self):

        if self._yasara_dir is None:
            raise Exception("yasara_dir has not been set")

        if self.execution_root_dir is None:
            raise Exception("execution root not set")

        if self.model_root_dir is None:
            raise Exception("model root not set")

        if self.template_blacklist is None:
            raise Exception("blacklist not set")

    def _add_template_to_blacklist (self, pdbid):

        open (self.template_blacklist, 'a').write ('%s\n' % pdbid.lower ())

    def _template_in_blacklist (self, pdbid):

        return (open (self.template_blacklist, 'r').read ().find (pdbid.lower ()) != -1)

    # Use this function to get the template sequences:
    def getChainOrderAndSeqs(self, tempobj):

        self._check_init()

        # chainOrder is the order that yasara maintains.

        chainOrder = self.yasara.ListMol('obj %i protein' % tempobj, 'MOL')
        s = self.yasara.SequenceMol('obj %i protein' % tempobj)

        uniqueChainOrder = []
        tempChainSeqs = {}
        for i in range(len(chainOrder)):
            chain = chainOrder[i]
            seq = ''

            for s in self.yasara.SequenceMol(
                    'obj %i mol %s protein' % (tempobj, chain)):
                seq += s

            # Some templates have two protein chains with the same ID:
            if chain not in tempChainSeqs:
                uniqueChainOrder.append (chain)

            tempChainSeqs[chain] = seq

        return [uniqueChainOrder, tempChainSeqs]

    def getChainSeq(self, tempobj, chain):
        self._check_init()

        return self.yasara.SequenceMol('%s obj %i protein' %
                                       (chain, tempobj))[0]

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
        self._check_init()

        self.yasara.Clear()

        tempobj = self.yasara.LoadPDB(tempac, download='best')[0]
        self.yasara.DelObj('not %i' % tempobj)

        # Count the number of molecules in unoligomerized state:
        nMolsUnoligomerized = \
            len(self.yasara.ListMol('obj %i protein' % tempobj, 'MOL'))

        # oligomerisation = 1
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
                tempobj = self.yasara.LoadPDB(tempac, download='best')[0]

            # Count the number of molecules in oligomerized state:
            nMolsOligomerized = len(
                self.yasara.ListMol('obj %i protein' % tempobj, 'MOL'))

            # We want as many chains as possible in the template.
            # If oligomerisation decreased the number of molecules, undo it!
            if nMolsOligomerized < nMolsUnoligomerized:

                self.yasara.Clear()
                tempobj = self.yasara.LoadPDB(tempac, download='best')[0]

        # If BuildSymRes throws an exception, we just continue
        try:
            self.yasara.BuildSymRes('obj %i' % tempobj)
        except Exception as e:
            _log.error('Cannot execute BuildSymRes on {}: {}'.format (tempac, e.args [0]))

        # Make sure there's only one chain for each chain identifier:
        chainOrder = self.yasara.ListMol('obj %i protein' % tempobj, 'MOL')
        for i in range (len (chainOrder)):
            if i > 0 and chainOrder [i - 1] == chainOrder [i]:
                self.yasara.JoinMol (
                    '%s and obj %i and Protein' %
                    (chainOrder[i], tempobj))
                

        self.yasara.CleanObj (tempobj)

        # Make sure there are no chains with sequence XXXXXXXXXXXXXXXXXX,
        # Otherwise, yasara would remove the entire chain.
        self.yasara.SwapRes (
            'Protein and UNK and atom CA and atom CB', 'ALA')
        self.yasara.SwapRes (
            'Protein and UNK and atom CA and not atom CB', 'GLY')

        # Delete buffer molecules:
        self.yasara.DelRes ('HOH H2O DOD D2O TIP WAT SOL ACT ACM ACY ' +
                            'EDO EOH FLC FMT GOL I3M IPA MLI MOH PEO ' +
                            'PO4 SO3 SO4 _TE UNX ACE')

        # Final checks for unfixable errors:
        chainOrder = self.yasara.ListMol ('obj %i protein' % tempobj, 'MOL')
        for chain in chainOrder:
            nocc = 0
            for c in chainOrder:
                if c == chain:
                    nocc += 1
            if nocc > 1:
                self._add_template_to_blacklist (tempac)
                raise Exception ("chain %s occurs more than once after cleaning" % chain)

        return [tempobj, nMolsOligomerized / nMolsUnoligomerized]


    # This function builds a model for one alignment.
    # It must be given a directory to work in,
    # the main target's id, sequence and domain range object,
    # a uniprot species id, a template ID (pdbid and chain)
    def _build_for_domain (self, modelDir, mainTargetID, uniprotSpeciesName,
                           mainTemplateID, mainTargetSeq, mainDomainRange):

            modelPath = os.path.join(modelDir, 'target.pdb')
            alignmentFastaPath = os.path.join(modelDir, 'align.fasta')
            selectedTargetsPath = \
                os.path.join(modelDir, 'selected-targets.txt')

            open(selectedTargetsPath, 'w') \
                .write('template: %s\n' % mainTemplateID.pdbac)

            # Load template and perform all necessary modifications:
            tempobj, oligomerisation = \
                self._set_template (mainTemplateID.pdbac)

            # One alignment per chain id:
            alignments = {}

            # Get the template's sequences:
            chainOrder, templateChainSequences = \
                self.getChainOrderAndSeqs(tempobj)

            # Remove funny insertions/deletions from main target:
            mainTargetSeq = adjustTargetSequence (mainTargetSeq, templateChainSequences [mainTemplateID.chainID], uniprotSpeciesName)

            _log.debug("template sequences after _set_template:\n " +
                       str(templateChainSequences))


            # A template can have more than one chain, but the main target is
            # not necessarily homologous to all of these chains. So we must make
            # sure that we find other targets for those other template chains.
            # However, we only need to do this for chains that actually interact
            # with the main target's template chain...

            time_start = time()

           
            # Part of the main target sequence that's in the domain range:
            mainDomainSeq = \
                mainTargetSeq [mainDomainRange.start: mainDomainRange.end]

            # Add the first alignments,
            # involving the main target sequence:
            # (the main target might be selected for more than one template sequence)
            for chainID in pickTemplateChainsFor (
                    templateChainSequences, mainTargetSeq):

                # Of coarse, we've already made an alignment previously
                # while determining the covered domain. However, we want
                # to make a new final alignment here, using a fine-tuned
                # program.

                # It's based on secondary structure, so we get that
                # information first:
                tempCAs, templateChainSeq, templateChainSecStr = \
                    getChainCAsSeqSecStr(self.yasara, tempobj, chainID)

                _log.debug ("aligning input target sequence \n{}"
                            .format(mainDomainSeq)
                            + "\nto template {} {}"
                            .format(mainTemplateID.pdbac, chainID)
                            + " sequence\n{}\n{}"
                            .format(templateChainSeq, templateChainSecStr))

                alignments [chainID] = aligner \
                    .msaAlign(templateChainSeq,
                              templateChainSecStr, mainDomainSeq)
                open (selectedTargetsPath, 'a') \
                    .write('\tmodeling main target %s on chain %s\n' %
                           (mainTargetID, chainID))

            # Try to find and align target sequences for
            # interacting chains in the template, while keeping
            # in mind which residues interact and must thus be
            # covered by the alignment.
            # We expand the set of involved template
            # chains with every iteration, until all template
            # chains have been added.
            yasaraChains = {}
            for c in self.yasara.ListMol ("obj %i and protein" % tempobj, "MOL"):
                yasaraChains [c] = YasaraChain (self.yasara, tempobj, c)
            while len (alignments) < \
                    len (self.yasara.ListMol (
                        'obj %i and protein' % tempobj)):

                # First, make python remember to which
                # chains the candidate chains interact:
                candidateChainInteractsWith = {}
                for c in alignments:
                    for chainID in \
                            listInteractingChains (yasaraChains [c]):

                        # Skip those, that we've already aligned
                        # to prevent infinite loops:
                        if chainID in alignments:
                            continue

                        if chainID not in candidateChainInteractsWith:
                            candidateChainInteractsWith[chainID] = []

                        candidateChainInteractsWith [chainID].append(c)

                if len (candidateChainInteractsWith) <= 0:

                    # No more interacting chains to add.
                    break

                # iterate over chains that might interact with the chains
                # that are already in the set:
                for chainID in candidateChainInteractsWith:

                    # Gather the alignments of
                    # the interaction partner chains.
                    # These are chains that we've already added to the set,
                    # so their alignments are already there.
                    interactingChainAlignments = {}
                    for interactingChainID in \
                            candidateChainInteractsWith[chainID]:

                        interactingChainAlignments [interactingChainID] = \
                            alignments [interactingChainID]

                    tempCAs, templateChainSeq, templateChainSecStr = \
                        getChainCAsSeqSecStr (self.yasara, tempobj, chainID)

                    # Find targets from the target species, homologous to this template chain:
                    potentialTargetSeqs = findTargets (templateChainSeq,
                                                       uniprotSpeciesName)

                    _log.debug ("found {} potential target sequences for {} {}"
                               .format (len (potentialTargetSeqs),
                                       mainTemplateID.pdbac, chainID))

                    yasaraChain = yasaraChains [chainID]

                    # Pick the target with the hightest sequence identity,
                    # and domain coverage.
                    targetsInterproRanges = {}
                    bestPID = 0.0
                    selectedTarget = None
                    for targetID in potentialTargetSeqs.keys():

                        _log.debug(
                            "aligning potential target sequence\n{}\n"
                            .format(potentialTargetSeqs[targetID])
                            + " to interacting template {} {}"
                            .format(mainTemplateID.pdbac, chainID)
                            + " sequence\n{}\n{}"
                            .format(templateChainSeq, templateChainSecStr))

                        # If we have targets with extremely high coverage,
                        # then we don't need interpro
                        alignment = aligner.msaAlign(
                            yasaraChain.seq,
                            yasaraChain.secstr,
                            potentialTargetSeqs[targetID])

                        pcov, pid = getCoverageIdentity(
                            alignment['template'], alignment['target'])

                        if pcov > 90.0: # Coverage is large enough, take the whole sequence:

                            if chainID not in alignments or pid > bestPID:
                                alignments[chainID] = alignment
                                bestPID = pid
                                selectedTarget = targetID

                        else: # no luck, bother interpro:
                            targetsInterproRanges [targetID] = interpro \
                                .getInterproDomainLocations (mainTargetSeq)

                    if chainID not in alignments: # Not been added yet at this point

                        _log.debug ("no 90 percent coverage target for {} {}"
                                   .format(mainTemplateID.pdbac, chainID))

                        # Find out which targets, given their domain ranges,
                        # preserve the interaction with the other chains when
                        # aligned to the template chain:
                        picker = InteractionPicker (chainID,
                            yasaraChains, interactingChainAlignments)
                        alignmentTriplesPerTarget = \
                            domainalign.pickAlignments(
                                yasaraChain, potentialTargetSeqs,
                                targetsInterproRanges, picker)

                        _log.debug("got {} targets passed by domainalign"
                                   .format(len(alignmentTriplesPerTarget))
                                   + " interacton picker for {} {}"
                                   .format(mainTemplateID.pdbac, chainID))

                        # Iterate over targets that passed:
                        for targetID in alignmentTriplesPerTarget:

                            nrange = len(alignmentTriplesPerTarget [targetID])
                            _log.debug("got {} interpro ranges for {} on {} {}"
                                       .format(nrange, targetID,
                                               mainTemplateID.pdbac, chainID))

                            # Join the alignments for this one target into one,
                            # maximizing the coverage of the template chain:
                            domain_alignment = domainalign \
                                .joinAlignmentsToBestTemplateCoverage(
                                    alignmentTriplesPerTarget [targetID])

                            # Remove gaps from aligned target sequence:
                            domain_target_seq = \
                                domain_alignment['target'].replace('-', '')

                            _log.debug(
                                "aligning {} domain target sequence\n{}"
                                .format(targetID, domain_target_seq)
                                + "\nto interacting template {} {}"
                                .format(mainTemplateID.pdbac, chainID)
                                + " sequence\n{}\n{}"
                                .format(templateChainSeq, templateChainSecStr))

                            # Realign once more using the complete target
                            # range that we found by joining the alignments:
                            alignment = aligner.msaAlign(
                                templateChainSeq, templateChainSecStr,
                                domain_target_seq)

                            # Determine the best target for this chain id:
                            # (highest identity)
                            nalign, pid = getNalignIdentity(
                                alignment ['target'], alignment ['template'])
                            if chainID not in alignments or pid > bestPID:
                                alignments [chainID] = alignment
                                bestPID = pid
                                selectedTarget = targetID

                    _log.debug("selected target for interacting {} {} is {}"
                               .format (mainTemplateID.pdbac, chainID,
                                       selectedTarget))

                    # Occasionally, we don't find any target for a template chain.
                    # It's solved here:
                    if chainID not in alignments:

                        _log.debug ("putting poly-A on chain {} of {}"
                                   .format (chainID, mainTemplateID.pdbac))

                        # Place alanines in the target sequence, with same
                        # length as templare chain:
                        alignments [chainID] = {
                            'template': yasaraChain.seq,
                            'target': 'A' * len(yasaraChain.seq)}
                        selectedTarget = 'poly-A'
                        self.yasara.OccupMol(chainID, 0.0)

                    # We document in a file which targets we used for this model:
                    open(selectedTargetsPath, 'a') \
                        .write('\tmodeling target %s on chain %s\n'
                               % (selectedTarget, chainID))

            time_targets = time()

            time_log ("took %d seconds to pick targets for template\n" % (time_targets - time_start))

            # < end of interaction finding iter

            # Delete chains that weren't aligned, assuming there's no
            # interaction with the main target's homologs:
            for chainID in chainOrder:
                if chainID not in alignments:

                    _log.debug("deleting not-interacting chain {} of {}"
                               .format(chainID, mainTemplateID.pdbac))

                    self.yasara.DelMol('%s and protein' % chainID)
                    open(selectedTargetsPath, 'a') \
                        .write('\tdeleting not-interacting '
                               + 'chain %s from template\n' % chainID)

            # The set of chains might have changed:
            chainOrder, templateChainSequences = \
                self.getChainOrderAndSeqs (tempobj)

            # Make the alignment file for yasara:
            writeAlignmentFasta(
                chainOrder, alignments, mainTemplateID.pdbac,
                alignmentFastaPath)

            # Start the modeling run:
            self.modelWithAlignment(alignmentFastaPath, tempobj)

            if not os.path.isfile ("target.yob"):
                raise Exception ("yasara did not build")

            # Save the model in PDB format:
            self.yasara.SavePDB (tempobj, modelPath)

            _log.info ("sucessfully created " + modelPath)

            time_model = time()

            time_log ("yasara modeling run took %d seconds\n" % (time_model - time_targets))


    # This is the main function for building a set of models.
    # Input is a target sequence and a uniprot species id.
    # The function seqrches its own templates. 
    #
    # Optional:
    #   requireRes: a residue number in target sequence (1,2,3,..) that must be covered by the model.
    #   overwrite: true if old models must be rebuilt, false if they must be skipped
    def modelProc (self, mainTargetSeq, uniprotSpeciesName, requireRes=None,
                  overwrite=False):
        self._check_init()

        if not os.path.isdir(self.model_root_dir):
            os.mkdir(self.model_root_dir)

        if not os.path.isdir(self.execution_root_dir):
            os.mkdir(self.execution_root_dir)

        # A model can have multiple target sequences, since the template
        # can have multiple chains.
        # The target sequence, that the model was built for, will be
        # called the 'main' target sequence.

        # Unique ID for this sequence:
        mainTargetID = idForSeq (mainTargetSeq)

        # Determine all domain ranges within our main target sequence.
        # (source: interpro)
        ranges = interpro.getInterproDomainLocations (mainTargetSeq)

        # yasara sticks to the directory where it was started,
        # so make a special directory for yasara to run in and
        # let it store all its output files there:
        runDir = os.path.join(self.execution_root_dir,
                              'run-yasara-%i' % (os.getpid()))
        if os.path.isdir(runDir):
            shutil.rmtree(runDir)
        os.mkdir(runDir)
        os.chdir(runDir)

        time_start = time()

        # Blast and filter alignments that satisfy the given set of
        # domain ranges. No halfway cut domains are alowed!
        mainTargetAlignments = \
            domainalign.getAlignments (ranges, mainTargetSeq)

        time_after_alignments = time()

        time_log ("took %i seconds to compute and filter alignments\n" % (time_after_alignments - time_start))

        if len (mainTargetAlignments) <= 0:
            _log.info ('no alignments found for sequence:\n' + mainTargetSeq)

        modelPaths = []
        failedModels = []

        realign = False

        # Iterate over all alignments that we've got. Any alignment
        # is a potential model.
        for mainDomainRange, mainTemplateID, mainDomainAlignment in \
                mainTargetAlignments:

            # skip ranges that don't cover the requireRes (if given)
            if requireRes and \
                    ((requireRes - 1) < mainDomainRange.start or
                     (requireRes - 1) >= mainDomainRange.end):
                continue

            modelname = '%s_%s_%i-%i' % \
                (mainTargetID, uniprotSpeciesName,
                 mainDomainRange.start + 1, mainDomainRange.end)

            modelDir = os.path.join (self.model_root_dir, modelname)
            modelArchive = modelDir + '.tgz'

            # A unique lockfile name, based on the model name:
            lockfile_path = modelDir + '_lock'

            # we dont want two threads building the same model, so lock:
            lock = FileLock (lockfile_path)

            with lock:
                if os.path.isfile (modelArchive) and not overwrite:

                # Found the archive for this model before starting, see if it's populated:

                    members = tarfile.open (modelArchive, mode='r:gz').getnames()
                    pdbpath = os.path.join (modelname, 'target.pdb')
                    if pdbpath in members:

                        modelPaths.append (modelArchive)

                        _log.info('%s exists, skipping..' % modelArchive)
                        continue

                # the model directory is going to be turned into an archive,
                # All created files will be sent there eventually,
                # but it's not yasara's working directory.
                if not os.path.isdir(modelDir):
                    os.mkdir(modelDir)

                modelPath = os.path.join(modelDir, 'target.pdb')

                if os.path.isfile (modelPath) and not overwrite:
                # Found the output directory for this model already before we started.
                # (just not yet archived)

                    # archive and clean up:
                    os.chdir(self.model_root_dir)
                    tf = tarfile.open(modelname + '.tgz', 'w:gz')
                    tf.add(modelname)  # refers to modelDir
                    tf.close()
                    shutil.rmtree(modelDir)

                    modelPaths.append(modelArchive)
                    _log.info('%s already exists, skipping..' % modelArchive)
                    continue

                # We're now sure that the model doesn't exist yet.
                # At this point we start the actual model building.
                try:
                    self._build_for_domain (modelDir, mainTargetID,
                                            uniprotSpeciesName, mainTemplateID,
                                            mainTargetSeq, mainDomainRange)

                except: # An error ocurred during modeling.

                    # Don't exit, just print error to log and
                    # move on to the next alignment.

                    exc_type, exc_value, exc_traceback = sys.exc_info()
                    stacktrace = ''.join(traceback.format_exception(
                        exc_type, exc_value, exc_traceback))

                    _log.error('an exception occured for {}:\n{}'
                               .format(modelname, stacktrace))


                    if self._template_in_blacklist (mainTemplateID.pdbac):

                        # The template has just been blacklisted, need to find a new one
                        realign = True
                    else:
                        failedModels.append (modelname)

                    # Also include error and yasara scene in the model dir,
                    # for debugging and yasara bug reports.
                    open (os.path.join (modelDir, 'errorexit.txt'), 'w') \
                        .write(stacktrace)
                    self.yasara.SaveSce (os.path.join (modelDir, 'errorexit.sce'))

                # At this point the model building has finished.
                # Move files from temporary run directory to their final destination..

                # Move all the files that yasara created:
                for f in os.listdir (runDir):
                    os.rename (os.path.join (runDir, f),
                               os.path.join (modelDir, f))

                # Create archive:
                parent = os.path.dirname (modelDir)
                os.chdir (parent)
                tf = tarfile.open (modelname + '.tgz', 'w:gz')
                tf.add (modelname)  # refers to modelDir
                tf.close ()
                shutil.rmtree (modelDir)
                modelPaths.append (modelArchive)
                os.chdir (runDir)

            # end of this iteration, move over to next range ...

        # Clean up all runtime files:
        os.chdir (self.execution_root_dir)
        shutil.rmtree (runDir)

        if realign:

            # Repeat once more with different templates:
            return self.modelProc (mainTargetSeq, uniprotSpeciesName,
                                   requireRes, overwrite)

        elif len (failedModels) > 0:

            raise Exception ("the following models have failed: " +
                             str(failedModels))
        else:
            return modelPaths


    # This function is called after yasara has finished builing the model.
    # It must perform the final checks to be sure that the model is not too bad.
    # If everything is OK, it must return nothing.
    # If something is wrong, it returns a string telling what's wrong.
    def checkmodel(self, modelobj, templateAC, expectedChains, targetAlignment,
                   templateAlignment):
        self._check_init()

        # Check that all chains are there:
        chainOrder, templateChainSequences = self.getChainOrderAndSeqs (modelobj)
        if len (chainOrder) < len (expectedChains):
            return "missing %i chains" % (len (expectedChains) - len (chainOrder))

        # Verify that the alignment is good enough:
        # include one extra iteration to finalize the last chain
        nalign = 0
        nid = 0
        nseq = 1
        # include one extra iteration to finalize the last chain
        for i in range (len (targetAlignment) + 1):

            if i >= len(targetAlignment) or targetAlignment[i] == '|':
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
                if targetAlignment[i].isalpha() and \
                        templateAlignment[i].isalpha():
                    nalign += 1
                    if targetAlignment[i].upper() == \
                            templateAlignment[i].upper():
                        nid += 1

        # Verify that packing quality didn't decrease to much,
        # relative to template:
        pdbname = 'model%i.pdb' % os.getpid()
        self.yasara.SavePDB(modelobj, pdbname)

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
            return 'whatif output for model: ' + str(qr)

        qr = whatifClient.service.PackingQualityMolecule(templateAC)
        if hasattr(qr[0], 'value'):
            templatePackingQuality = qr[0].value
            open('packing.log', 'a') \
                .write('template packing quality is %.3f\n' %
                       templatePackingQuality)
        else:
            return 'whatif output for template: ' + str(qr)

        if modelPackingQuality < (templatePackingQuality - 1.0):
            return ("Model packing quality %.3f too low" +
                    ", compared to %.3f from template") % \
                (modelPackingQuality, templatePackingQuality)

    # This function starts the model building process by calling
    # yasara. 
    def modelWithAlignment (self, alignmentFastaPath, tempobj):
        self._check_init()

        self.yasara.Processors(1)

        # self.yasara.Clear()
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

modeler = Modeler()


def selectMostIdentical(targetSeq, seqs):
    best = None

    naligns = {}
    pids = {}
    for key in seqs.keys():

        d = {key: seqs[key], 'target': targetSeq}
        aligned = aligner.clustalAlign(d)

        nid = 0
        naligns[key] = 0
        for i in range(len(aligned['target'])):
            taa = aligned['target'][i]
            saa = aligned[key][i]

            if taa.isalpha() and saa.isalpha():
                naligns[key] += 1
                if taa.upper() == saa.upper():
                    nid += 1

        pids[key] = nid / naligns[key]

        if not best or pids[key] > pids[best]:
            best = key

    return best


def selectLongestAligned(hits, minLen=0):
    bestHitID = None
    bestHitAlignment = None

    for hitID in hits.keys():
        for alignment in hits[hitID]:
            identity = alignment.getIdentity()
            nalign = alignment.getNumberResiduesAligned()

            if nalign < minLen:
                continue

            if not bestHitID or \
                    nalign > bestHitAlignment.getNumberResiduesAligned() or \
                    nalign == bestHitAlignment.getNumberResiduesAligned() and \
                    identity > bestHitAlignment.getIdentity():
                bestHitID = hitID
                bestHitAlignment = alignment

    return [bestHitID, bestHitAlignment]


# This function finds uniprot isoforms for a given sequence.
# (a sequence that looks almost like it, except for some deletions)
# It makes sure the isoform comes from the given species.
def getIsoforms (seq, uniprotSpeciesName):
    hits = blaster.speciesBlast(seq, uniprotSpeciesName)

    selected = {}
    for hitID in hits:
        for ali in hits[hitID]:
            ac = hitID.split('|')[1]
            if ali.getIdentity() >= 99.0:
                if ac not in selected or \
                        selected[ac].getNumberResiduesAligned() < \
                        ali.getNumberResiduesAligned():
                    selected[ac] = ali

    return selected

# This function fills in deletions with pieces from isoforms in uniprot.
# It's necessary sometimes, because deletions might destroy the
# hydrofobic core of a model.
def adjustTargetSequence (targetseq, templateseq, uniprotSpeciesName):

    seqs = {}
    tarkey = 'target'
    seqs[tarkey] = targetseq
    tempkey = 'template'
    seqs[tempkey] = templateseq

    isoformkeys = []
    isoformHits = getIsoforms(targetseq, uniprotSpeciesName)
    for ac in isoformHits.keys():
        hit = isoformHits[ac]
        seqs[ac] = hit.subjectalignment.replace ('-', '').upper()
        if len(seqs[ac]) < 10000:  # avoid titin
            isoformkeys.append(ac)

    aligned = aligner.clustalAlign({
        tempkey: seqs[tempkey], tarkey: seqs[tarkey]})
    deletionRanges = identifyDeletedRegions(aligned[tarkey])

    if len (deletionRanges) <= 0:
        return targetseq

    adjusted = aligned[tarkey]

    # insertions must be added in reverse order, so that indices don't get
    # messed up
    deletionRanges.reverse()
    for i, f in deletionRanges:

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
                aligner.clustalAlign({tarkey: seqs[tarkey], key: seqs[key]})

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
                ali = aligner.clustalAlign({
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

    return adjusted.replace('-', '')


# For a given sequence, finds an ortholog in the requested species:
def findOrthologsToSeq (seq, uniprotSpeciesName):
    orthologs = {}

    hits = blaster.speciesBlast(seq, uniprotSpeciesName)
    for hitID in hits.keys():
        for alignment in hits[hitID]:

            pcov, pid = getCoverageIdentity(
                alignment.queryalignment, alignment.subjectalignment)

            if pid > 70.0:

                ac = hitID.split('|')[1]
                oseq = getUniprotSeq(ac)

                if pcov > 90.0:
                    orthologs[ac] = oseq

    _log.debug("found {} {} orthologs to\n{}"
               .format(len(orthologs), uniprotSpeciesName, seq))

    return orthologs


def findTargets(tempChainSeq, uniprotSpeciesName):
    return findOrthologsToSeq(tempChainSeq, uniprotSpeciesName)


def groupIdenticals(d):  # arg is dictionary of sequences
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
            pair = pair[0]+'_'+pair[1]

            # Aligning them all in one run can be very time-consuming

            if pair not in alignments:
                alignments[pair] = aligner.clustalAlign({
                    ID: d[ID],
                    otherID: d[otherID]
                })
            alignedD = alignments[pair]

            if getPercentageIdentity(alignedD[ID], alignedD[otherID]) >= 99.0:
                grouped[-1].append(otherID)
                ids.remove(otherID)

    return grouped


def pickTemplateChainsFor(chainSequences, targetSeq):
    identicalChains = groupIdenticals(chainSequences)

    hitChain = selectHighestSeqID(targetSeq,
                                  filterMinLength(30, chainSequences))

    for group in identicalChains:

        if hitChain in group:
            return group

    raise Exception('no hit found for %s in template' % targetSeq)
