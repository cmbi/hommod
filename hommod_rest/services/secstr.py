import os
import sys
from hommod_rest.services.modelutils import parseDSSP

# Secondary structure can either be taken from dssp or yasara.
# Use yasara when dssp is not available.
class SecStrProvider(object):

    def __init__(self):

        self.dssp_dir = ''
        self._yasara_dir = None

    def _check_config(self):

        if not self.dssp_dir:
            raise Exception("dssp dir not set")
        if not self._yasara_dir:
            raise Exception("yasara dir not set")

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

    def hasSecStr(self, template):

        self._check_config()

        dssppath = '%s/%s.dssp' % (self.dssp_dir, template.pdbac.lower())
        if os.path.isfile(dssppath):

            d = parseDSSP(dssppath)
            return template.chainID in d
        else:
            obj = self.yasara.LoadPDB(template.pdbac, download='Yes')[0]
            for ss in self.yasara.SecStrRes('obj %i' % obj):
                if ss != 'X':
                    self.yasara.DelObj(obj)
                    return True
            self.yasara.DelObj(obj)


secstr = SecStrProvider()
