import logging
import os
import re
import shutil
import sys
import tarfile
import traceback

from suds import client

from modelutils import (getChainCAsSeqSecStr, getNalignIdentity,
                        getCoverageIdentity, YasaraChain,
                        identifyDeletedRegions, getPercentageIdentity,
                        filterMinLength, minIdentity,
                        getUniprotSeq, idForSeq)

from interaction import listInteractingChains, InteractionPicker

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


def selectHighestSeqID(seq, fromd):
    qid = 'query'

    bestID = 0
    best = None

    for homolog in fromd.keys():

        # Aligning them all in one run can be very time consuming !

        aligned = aligner.clustalAlign({qid: seq, homolog: fromd[homolog]})
        pid = getPercentageIdentity(aligned[qid], aligned[homolog])

        if pid > bestID:
            bestID = pid
            best = homolog

    return best


def onlyX(seq):
    return (len(seq) == seq.count('X'))


def withoutHomomerDuplicates(chainOrder, chainSequences):
    count = {}
    seq2chains = {}

    for chain in chainOrder:
        seq = chainSequences[chain]
        if seq not in count:
            count[seq] = 0
            seq2chains[seq] = []
        count[seq] += 1
        seq2chains[seq].append(chain)

    # Check that every sequence occurs the same number of times:
    sameNum = True
    for n in count.values()[1:]:
        if n != count.values()[0]:
            sameNum = False
            break

    if sameNum:
        newChainSequences = {}
        newChainOrder = []
        for seq in count.keys():
            newChainSequences[seq2chains[seq][0]] = seq

        for chain in chainOrder:
            if chain in newChainSequences:
                newChainOrder.append(chain)

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
    def __init__(self, yasara_dir=None):
        self._yasara_dir = yasara_dir
        self.execution_root_dir = None
        self.model_root_dir = None

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
        self.yasara.info.mode = 'txt'

    def _check_init(self):
        if self._yasara_dir is None:
            raise Exception("yasara_dir has not been set")

        if self.execution_root_dir is None:
            raise Exception("execution root not set")

        if self.model_root_dir is None:
            raise Exception("model root not set")

    def getChainOrderAndSeqs(self, tempobj):

        self._check_init()

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
                uniqueChainOrder.append(chain)

            tempChainSeqs[chain] = seq

        return [uniqueChainOrder, tempChainSeqs]

    def getChainSeq(self, tempobj, chain):
        self._check_init()

        return self.yasara.SequenceMol('%s obj %i protein' %
                                       (chain, tempobj))[0]

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
            _log.error('Cannot execute BuildSymRes on {}: {}' % (tempac, e))

        # Make sure there's only one chain for each chain identifier:
        chainOrder = self.yasara.ListMol('obj %i protein' % tempobj, 'MOL')
        for i in range(len(chainOrder)):
            if i > 0 and chainOrder[i - 1] == chainOrder[i]:
                residues = self.yasara.ListRes(
                    'mol %s Protein' % chainOrder[i],
                    "RESNUM")
                self.yasara.JoinMol(
                    '%s res %s obj %i Protein' %
                    (chainOrder[i], residues[-1], tempobj))

        self.yasara.CleanObj(tempobj)

        # Make sure there are no chains with sequence XXXXXXXXXXXXXXXXXX,
        # Otherwise, yasara would remove the entire chain.
        self.yasara.SwapRes(
            'Protein and UNK and atom CA and atom CB', 'ALA')
        self.yasara.SwapRes(
            'Protein and UNK and atom CA and not atom CB', 'GLY')

        self.yasara.DelRes('HOH H2O DOD D2O TIP WAT SOL ACT ACM ACY ' +
                           'EDO EOH FLC FMT GOL I3M IPA MLI MOH PEO ' +
                           'PO4 SO3 SO4 _TE UNX ACE')

        return [tempobj, nMolsOligomerized / nMolsUnoligomerized]

    def separateChains(self, obj, chainsToSeparate):
        self._check_init()

        chains = self.yasara.ListMol("obj %i and protein" % obj)
        chainsToKeep = []
        for chain in chains:
            if chain not in chainsToSeparate:
                chainsToKeep.append(chain)

        newobj = self.yasara.DuplicateObj(obj)[0]

        for chain in chainsToKeep:
            self.yasara.DelMol('obj %i and %s' % (newobj, chain))

        for chain in chainsToSeparate:
            self.yasara.DelMol('obj %i and %s' % (obj, chain))

        return newobj

    def _build_for_domain(self, modelDir, mainTargetID, uniprotSpeciesName,
                          mainTemplateID, mainTargetSeq, mainDomainRange):

            modelPath = os.path.join(modelDir, 'target.pdb')
            alignmentFastaPath = os.path.join(modelDir, 'align.fasta')
            selectedTargetsPath = \
                os.path.join(modelDir, 'selected-targets.txt')

            open(selectedTargetsPath, 'w') \
                .write('template: %s\n' % mainTemplateID.pdbac)

            # Load template:
            tempobj, oligomerisation = \
                self._set_template(mainTemplateID.pdbac)

            alignments = {}

            # Get the template's sequences:
            chainOrder, templateChainSequences = \
                self.getChainOrderAndSeqs(tempobj)

            _log.debug("template sequences after _set_template:\n " +
                       str(templateChainSequences))

            # Add the first alignments,
            # belonging to the input target protein
            mainDomainSeq = \
                mainTargetSeq[mainDomainRange.start: mainDomainRange.end]
            for chainID in pickTemplateChainsFor(
                    templateChainSequences, mainTargetSeq):

                tempCAs, templateChainSeq, templateChainSecStr = \
                    getChainCAsSeqSecStr(self.yasara, tempobj, chainID)

                _log.debug("aligning input target sequence \n{}"
                           .format(mainDomainSeq)
                           + "\nto template {} {}"
                           .format(mainTemplateID.pdbac, chainID)
                           + " sequence\n{}\n{}"
                           .format(templateChainSeq, templateChainSecStr))

                alignments[chainID] = aligner \
                    .msaAlign(templateChainSeq,
                              templateChainSecStr, mainDomainSeq)
                open(selectedTargetsPath, 'a') \
                    .write('\tmodeling main target %s on chain %s\n' %
                           (mainTargetID, chainID))

            # Try to find and align targets for
            # interacting chains in the template:
            while len(alignments) < \
                    len(self.yasara.ListMol(
                        'obj %i and protein' % tempobj)):

                # Make python remember to which
                # chains the candidate chains interact
                candidateChainInteractsWith = {}
                for c in alignments.keys():
                    for chainID in \
                            listInteractingChains(
                                YasaraChain(self.yasara, tempobj, c)):

                        # Skip those, that we've already aligned
                        # to prevent infinite loops:
                        if chainID in alignments:
                            continue

                        if chainID not in candidateChainInteractsWith:
                            candidateChainInteractsWith[chainID] = []

                        candidateChainInteractsWith[chainID].append(c)

                # iterate over chains that might interact with the chains
                # already aligned
                for chainID in candidateChainInteractsWith.keys():

                    # Gather the alignments of
                    # the interaction partner chains
                    interactingChainAlignments = {}
                    for interactingChainID in \
                            candidateChainInteractsWith[chainID]:
                        interactingChainAlignments[interactingChainID] = \
                            alignments[interactingChainID]

                    tempCAs, templateChainSeq, templateChainSecStr = \
                        getChainCAsSeqSecStr(self.yasara, tempobj, chainID)

                    potentialTargetSeqs = findTargets(templateChainSeq,
                                                      uniprotSpeciesName)

                    _log.debug("found {} potential target sequences for {} {}"
                               .format(len(potentialTargetSeqs),
                                       mainTemplateID.pdbac, chainID))

                    yasaraChain = YasaraChain(self.yasara, tempobj, chainID)

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
                        if pcov > 90.0:

                            if chainID not in alignments or pid > bestPID:
                                alignments[chainID] = alignment
                                bestPID = pid
                                selectedTarget = targetID

                        else:
                            targetsInterproRanges[targetID] = interpro \
                                .getInterproDomainLocations(mainTargetSeq)

                    if chainID not in alignments:

                        _log.debug("no 90 percent coverage target for {} {}"
                                   .format(mainTemplateID.pdbac, chainID))

                        picker = InteractionPicker(
                            yasaraChain, interactingChainAlignments)
                        alignmentTriplesPerTarget = \
                            domainalign.pickAlignments(
                                yasaraChain, potentialTargetSeqs,
                                targetsInterproRanges, picker)

                        _log.debug("got {} targets passed by domainalign"
                                   .format(len(alignmentTriplesPerTarget))
                                   + " interacton picker for {} {}"
                                   .format(mainTemplateID.pdbac, chainID))

                        for targetID in alignmentTriplesPerTarget.keys():

                            nrange = len(alignmentTriplesPerTarget[targetID])
                            _log.debug("got {} interpro ranges for {} on {} {}"
                                       .format(nrange, targetID,
                                               mainTemplateID.pdbac, chainID))

                            domain_alignment = domainalign \
                                .joinAlignmentsToBestTemplateCoverage(
                                    alignmentTriplesPerTarget[targetID])

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
                            # range that we found:
                            alignment = aligner.msaAlign(
                                templateChainSeq, templateChainSecStr,
                                domain_target_seq)

                            nalign, pid = getNalignIdentity(
                                alignment['target'], alignment['template'])
                            if chainID not in alignments or pid > bestPID:
                                alignments[chainID] = alignment
                                bestPID = pid
                                selectedTarget = targetID

                    _log.debug("selected target for interacting {} {} is {}"
                               .format(mainTemplateID.pdbac, chainID,
                                       selectedTarget))

                    if chainID not in alignments:

                        _log.debug("putting poly-A on chain {} of {}"
                                   .format(chainID, mainTemplateID.pdbac))

                        alignments[chainID] = {
                            'template': yasaraChain.seq,
                            'target': 'A' * len(yasaraChain.seq)}
                        selectedTarget = 'poly-A'
                        self.yasara.OccupMol(chainID, 0.0)

                    open(selectedTargetsPath, 'a') \
                        .write('\tmodeling target %s on chain %s\n'
                               % (selectedTarget, chainID))

            # < end of interaction finding iter

            # Delete chains that weren't aligned, assuming there's no
            # interaction
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
                self.getChainOrderAndSeqs(tempobj)

            # Make the alignment file for yasara:
            writeAlignmentFasta(
                chainOrder, alignments, mainTemplateID.pdbac,
                alignmentFastaPath)

            # Start the modeling run:
            self.modelWithAlignment(alignmentFastaPath, tempobj)

            # Save the model in PDB format:
            self.yasara.SavePDB(tempobj, modelPath)

            _log.info("sucessfully created " + modelPath)

    def modelProc(self, mainTargetSeq, uniprotSpeciesName, requireRes=None,
                  overwrite=False):
        self._check_init()

        if not os.path.isdir(self.model_root_dir):
            os.mkdir(self.model_root_dir)

        if not os.path.isdir(self.execution_root_dir):
            os.mkdir(self.execution_root_dir)

        mainTargetID = idForSeq(mainTargetSeq)

        ranges = interpro.getInterproDomainLocations(mainTargetSeq)

        # yasara sticks to the directory where it was started,
        # so make a special directory for yasara to run in and
        # let it store all its output files there:
        runDir = os.path.join(self.execution_root_dir,
                              'run-yasara-%i' % (os.getpid()))
        if os.path.isdir(runDir):
            shutil.rmtree(runDir)
        os.mkdir(runDir)
        os.chdir(runDir)

        mainTargetAlignments = \
            domainalign.getAlignments(ranges, mainTargetSeq)

        if len(mainTargetAlignments) <= 0:
            _log.info('no alignments found for sequence:\n' + mainTargetSeq)

        modelPaths = []
        failedModels = []

        for mainDomainRange, mainTemplateID, mainDomainAlignment in \
                mainTargetAlignments:
            # skip ranges that don't cover
            if requireRes and \
                    ((requireRes - 1) < mainDomainRange.start or
                     (requireRes - 1) >= mainDomainRange.end):
                continue

            modelname = '%s_%s_%i-%i' % \
                (mainTargetID, uniprotSpeciesName,
                 mainDomainRange.start + 1, mainDomainRange.end)

            modelDir = os.path.join(self.model_root_dir, modelname)
            modelArchive = modelDir + '.tgz'

            if os.path.isfile(modelArchive) and not overwrite:

                members = tarfile.open(modelArchive, mode='r:gz').getnames()
                pdbpath = os.path.join(modelname, 'target.pdb')
                if pdbpath in members:

                    modelPaths.append(modelArchive)

                    _log.info('%s exists, skipping..' % modelArchive)
                    continue

            # the model directory is going to be turned into an archive,
            # All created files will be sent there eventually,
            # but it's not yasara's working directory.
            if not os.path.isdir(modelDir):
                os.mkdir(modelDir)

            modelPath = os.path.join(modelDir, 'target.pdb')

            if os.path.isfile(modelPath) and not overwrite:

                # archive and clean up:
                os.chdir(self.model_root_dir)
                tf = tarfile.open(modelname + '.tgz', 'w:gz')
                tf.add(modelname)  # refers to modelDir
                tf.close()
                shutil.rmtree(modelDir)

                modelPaths.append(modelArchive)
                _log.info('%s already exists, skipping..' % modelArchive)
                continue

            try:
                self._build_for_domain(modelDir, mainTargetID,
                                       uniprotSpeciesName, mainTemplateID,
                                       mainTargetSeq, mainDomainRange)
            except:
                exc_type, exc_value, exc_traceback = sys.exc_info()
                stacktrace = ''.join(traceback.format_exception(
                    exc_type, exc_value, exc_traceback))

                _log.error('an exception occured for {}:\n{}'
                           .format(modelname, stacktrace))
                failedModels.append(modelname)

                open(os.path.join(modelDir, 'errorexit.txt'), 'w') \
                    .write(stacktrace)
                self.yasara.SaveSce(os.path.join(modelDir, 'errorexit.sce'))

            # Move all the files that yasara created:
            for f in os.listdir(runDir):
                os.rename(os.path.join(runDir, f), os.path.join(modelDir, f))

            parent = os.path.dirname(modelDir)
            os.chdir(parent)
            tf = tarfile.open(modelname + '.tgz', 'w:gz')
            tf.add(modelname)  # refers to modelDir
            tf.close()
            shutil.rmtree(modelDir)
            modelPaths.append(modelArchive)
            os.chdir(runDir)

            # end of this iteration, move over to next range ...

        os.chdir(self.execution_root_dir)
        shutil.rmtree(runDir)

        if len(failedModels) > 0:

            raise Exception("the following models have failed: " +
                            str(failedModels))
        else:
            return modelPaths

    def checkmodel(self, modelobj, templateAC, expectedChains, targetAlignment,
                   templateAlignment):
        self._check_init()

        # Check that all chains are there:
        chainOrder, templateChainSequences = self.getChainOrderAndSeqs(modelobj)
        if len(chainOrder) < len(expectedChains):
            return "missing %i chains" % (len(expectedChains) - len(chainOrder))

        # Verify that the alignment is good enough:
        # include one extra iteration to finalize the last chain
        nalign = 0
        nid = 0
        nseq = 1
        # include one extra iteration to finalize the last chain
        for i in range(len(targetAlignment) + 1):

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

    def modelWithAlignment(self, alignmentFastaPath, tempobj):
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
            speed='slow',
            loopsamples=20,
            resultfile='target'
        )
        self.yasara.Experiment("On")
        self.yasara.Wait("Expend")

    def modelWithSequence(self, targetFastaName, tempobj):
        self._check_init()

        self.yasara.Processors(1)

        self.yasara.Clear()
        self.yasara.ExperimentHomologyModeling(
            sequencefile=targetFastaName,
            templates="1, sameseq = 1",
            alignments=1,
            termextension=0,
            oligostate=32,
            looplenmax=10,
            animation='fast',
            speed='slow',
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
        aligned = aligner.alignClustal(d)

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


def getIsoforms(seq, uniprotSpeciesName):
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


def adjustTargetSequence(targetseq, templateseq, uniprotSpeciesName):
    # This function fills in deletions with pieces from isoforms in uniprot.
    seqs = {}
    tarkey = 'target'
    seqs[tarkey] = targetseq
    tempkey = 'template'
    seqs[tempkey] = templateseq

    isoformkeys = []
    isoformHits = getIsoforms(targetseq, uniprotSpeciesName)
    for ac in isoformHits.keys():
        hit = isoformHits[ac]
        seqs[ac] = hit.subjectalignment.replace('-', '').upper()
        if len(seqs[ac]) < 10000:  # avoid titin
            isoformkeys.append(ac)

    aligned = aligner.alignClustal({
        tempkey: seqs[tempkey], tarkey: seqs[tarkey]})
    deletionRanges = identifyDeletedRegions(aligned[tarkey])

    if len(deletionRanges) <= 0:
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
                aligner.alignClustal({tarkey: seqs[tarkey], key: seqs[key]})

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
                ali = aligner.alignClustal({
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


def findOrthologsToSeq(seq, uniprotSpeciesName):
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
