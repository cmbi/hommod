import sys
import imp
import platform
import os


class YasaraContext:
    def __init__(self, yasara_dir):
        sys.path.append(os.path.join(yasara_dir, 'pym'))
        sys.path.append(os.path.join(yasara_dir, 'plg'))
        self._yasara_module = imp.load_module('yasaramodule', *imp.find_module('yasaramodule'))

        self._start_yasara(yasara_dir)

    def _start_yasara(self, yasara_dir):
        self._com = self._yasara_module.yasara_communicator()

        executable = os.path.join(yasara_dir, "yasara")
        arglist = [executable, "-pym", str(self._com.port), "-txt"]

        self._pid = os.spawnv(os.P_NOWAIT, executable, arglist)
        self._com.accept()

    def _execute(self, messagedata):
        self._com.sendmessage(self._com.EXECUTE, messagedata)
        result = self._com.receivemessage(self._com.RESULT)
        return result

    def _run(self, command):
        if (command.find("\n") == -1):
            return self._execute(command)

        # Multiline command, send each line separately to catch all errors.
        return [self._execute(c) for c in command.split("\n")]

    def __delete__(self, instance):
        self._execute("Exit")
        os.waitpid(self._pid, 0)

    def ListRes(self, selection, format_=None):
        command = 'ListRes ' + self._yasara_module.selstr(selection) + ','
        if format_ is not None:
            command += 'Format=' + self._yasara_module.cstr(format_) + ','
        return self._run(command[:-1])

    def ListAtom(self, selection, format_=None):
        command = 'ListAtom ' + self._yasara_module.selstr(selection) + ','
        if format_ is not None:
            command += 'Format=' + self._yasara_module.cstr(format_) + ','
        return self._run(command[:-1])

    def ListMol(self, selection, format_=None):
        command = 'ListMol ' + self._yasara_module.selstr(selection) + ','
        if format_ is not None:
            command += 'Format=' + self._yasara_module.cstr(format_) + ','
        return self._run(command[:-1])

    def DownloadPDB(self, pdbid):
        command = 'LoadPDB Filename=%s, Download=yes' % self._yasara_module.cstr(pdbid)
        return self._run(command)

    def CD(self, dir_path):
        command = 'CD %s' % self._yasara_module.cstr(dir_path)
        return self._run(command)

    def SecStrRes(self, selection):
        command = 'SecStrRes %s' % self._yasara_module.selstr(selection)
        return self._run(command)

    def OligomerizeObj(self, selection):
        command = 'OligomerizeObj %s' % self._yasara_module.selstr(selection)
        return self._run(command)

    def BuildSymRes(self, selection):
        command = 'BuildSymRes %s' % self._yasara_module.selstr(selection)
        return self._run(command)

    def ExperimentHomologyModeling(self, alignfile=None,
                                         templateobj=None,
                                         templates=None,
                                         alignments=None,
                                         termextension=None,
                                         oligostate=None,
                                         looplenmax=None,
                                         animation=None,
                                         speed=None,
                                         loopsamples=None,
                                         resultfile=None,
                                         psiblasts=None):
        command = 'Experiment HomologyModeling\n'
        if alignfile is not None:
            command += '  alignfile %s\n' % self._yasara_module.cstr(alignfile)
        if templateobj is not None:
            command += '  templateobj %s\n' % self._yasara_module.cstr(templateobj)
        if templates is not None:
            command += '  templates %s\n' % self._yasara_module.cstr(templates)
        if alignments is not None:
            command += '  alignments %s\n' % self._yasara_module.cstr(alignments)
        if termextension is not None:
            command += '  termextension %s\n' % self._yasara_module.cstr(termextension)
        if oligostate is not None:
            command += '  oligostate %s\n' % self._yasara_module.cstr(oligostate)
        if looplenmax is not None:
            command += '  looplenmax %s\n' % self._yasara_module.cstr(looplenmax)
        if animation is not None:
            command += '  animation %s\n' % self._yasara_module.cstr(animation)
        if speed is not None:
            command += '  speed %s\n' % self._yasara_module.cstr(speed)
        if loopsamples is not None:
            command += '  loopsamples %s\n' % self._yasara_module.cstr(loopsamples)
        if resultfile is not None:
            command += '  resultfile %s\n' % self._yasara_module.cstr(resultfile)
        if psiblasts is not None:
            command += '  psiblasts %s\n' % self._yasara_module.cstr(psiblasts)

        self._run(command[:-1])

    def Experiment(self, name):
        self._run("Experiment %s" % name)

    def Wait(self, steps, unit=None):
        command='Wait %s,' % steps
        if (unit is not None):
            command += 'Unit=%s,' % cstr(unit)

        result = self._run(command[:-1])
        if result is not None and len(result) > 0:
            return result[0]
        return None

    def Processors(self, cputhreads=None, gpu=None):
        command = 'Processors '
        if cputhreads is not None:
            command += 'CPUThreads=%s,' % self._yasara_module.cstr(cputhreads)
        if gpu is not None:
            command += 'GPU=%s,' % self._yasara_module.cstr(gpu)

        return self._run(command[:-1])

    def Clear(self):
        self._run('Clear')

    def DelObj(self, selection):
        command = 'DelObj %s' % self._yasara_module.selstr(selection)
        self._run(command)

    def DelMol(self, selection):
        command = 'DelMol %s' % self._yasara_module.selstr(selection)
        self._run(command)

    def DelRes(self, selection):
        command = 'DelRes %s' % self._yasara_module.selstr(selection)
        self._run(command)

    def DelAtom(self, selection):
        command = 'DelAtom %s' % self._yasara_module.selstr(selection)
        self._run(command)

    def JoinMol(self, selection):
        command = 'JoinMol %s' % self._yasara_module.selstr(selection)
        self._run(command)

    def AddBond(self, selection1, selection2, order=None, update=None, lenmax=None):
        command = 'AddBond %s, %s,' % (self._yasara_module.selstr(selection1),
                                       self._yasara_module.selstr(selection2))
        if (order is not None):
            command += 'Order=%s,' % self._yasara_module.cstr(order)
        if (update is not None):
            command += 'Update=%s,' % self._yasara_module.cstr(update)
        if (lenmax is not None):
            command += 'LenMax=%s,' % self._yasara_module.cstr(lenmax)

        self._run(command[:-1])

    def SwapRes(self, selection, new, isomer=None):
        command = 'SwapRes %s,new=%s,' % (self._yasara_module.selstr(selection),
                                          self._yasara_module.cstr(new))
        if isomer is not None:
            command+='Isomer=%s,' % self._yasara_module.cstr(isomer)

        self._run(command)

    def CleanObj(self, selection):
        command = 'CleanObj %s' % self._yasara_module.selstr(selection)
        self._run(command)

    def SavePDB(self, selection, filename):
        command = 'SavePDB %s, Filename=%s' % (self._yasara_module.selstr(selection),
                                               self._yasara_module.cstr(filename))
        self._run(command)
