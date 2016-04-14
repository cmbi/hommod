import os
import sys
from hommod_rest.services.modelutils import parseDSSP, downloadPDB

import logging
_log = logging.getLogger(__name__)

# Secondary structure can either be taken from dssp or yasara.
# Use yasara when dssp is not available.
class SecStrProvider(object):

    def __init__(self):

        self.dssp_dir = ''
        self._yasara_dir = None

    def _check_config(self):

        from flask import current_app as flask_app

        self.dssp_dir = flask_app.config ['DSSPDIR']
        self.yasara_dir = flask_app.config ['YASARADIR']

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

        _log.info ("checking if %s_%s has secondary structure" % (template.pdbac, template.chainID))

        self._check_config()

        dssppath = '%s/%s.dssp' % (self.dssp_dir, template.pdbac.lower())
        if os.path.isfile(dssppath):

            d = parseDSSP(dssppath)

            _log.info ("from %s:\n%s" % (dssppath, str (d)))

            return template.chainID in d
        else:
            pdbfile = "%s.pdb" % template.pdbac
            open (pdbfile, 'w').write (downloadPDB (template.pdbac))
            obj = self.yasara.LoadPDB (pdbfile)[0]
            for ss in self.yasara.SecStrRes('obj %i' % obj):
                if ss != 'X':
                    self.yasara.DelObj(obj)

                    _log.info ("yasara reported secondary structure for %s_%s" % (template.pdbac, template.chainID))

                    return True
            self.yasara.DelObj(obj)

            _log.info ("yasara reported no secondary structure for %s_%s" % (template.pdbac, template.chainID))

            return False

secstr = SecStrProvider()
