import os
import sys
from hommod_rest.services.modelutils import parseDSSP, downloadPDB

import logging
_log = logging.getLogger(__name__)


class SecondaryStructureProvider(object):
    """
    Secondary structure can either be taken from dssp or yasara.
    Use yasara when dssp is not available.
    """
    def __init__(self, dssp_dir=None, yasara_dir=None):

        self._dssp_dir = dssp_dir
        self._yasara_dir = yasara_dir

    @property
    def dssp_dir(self):
        return self._dssp_dir

    @dssp_dir.setter
    def dssp_dir(self, dssp_dir):
        self._dssp_dir = dssp_dir

    @property
    def yasara_dir(self):
        return self._yasara_dir

    @yasara_dir.setter
    def yasara_dir(self, yasara_dir):
        self._yasara_dir = yasara_dir

        if os.path.isdir(yasara_dir):
            sys.path.append(os.path.join(yasara_dir, 'pym'))
            sys.path.append(os.path.join(yasara_dir, 'plg'))

            import yasaramodule
            self.yasara = yasaramodule
            self.yasara.info.mode = 'txt'

    def _check_config(self):
        if not self._dssp_dir:
            raise Exception("dssp_dir is not set")

        if not self._yasara_dir:
            raise Exception("yasara_dir is not set")

    def has_secondary_structure(self, template):
        _log.info("checking if %s_%s has secondary structure" % (template.pdbac, template.chainID))

        self._check_config()

        dssppath = '%s/%s.dssp' % (self.dssp_dir, template.pdbac.lower())
        if os.path.isfile(dssppath):

            d = parseDSSP(dssppath)

            _log.info("from %s:\n%s" % (dssppath, str(d)))

            return template.chainID in d
        else:
            obj = self.yasara.LoadPDB(template.pdbac, download='yes')[0]
            for ss in self.yasara.SecStrRes('obj %i and mol %s'
                                            % (obj, template.chainID)):
                if ss != 'X':
                    self.yasara.DelObj(obj)

                    _log.info("yasara reported secondary structure for %s_%s" % (template.pdbac, template.chainID))

                    return True
            self.yasara.DelObj(obj)

            _log.info("yasara reported no secondary structure for %s_%s" % (template.pdbac, template.chainID))

            return False

    def get_sequence_secondary_structure(self, template):
        self._check_config()

        dssppath = '%s/%s.dssp' % (self.dssp_dir, template.pdbac.lower())
        if os.path.isfile(dssppath):
            d = parseDSSP(dssppath)
            if template.chainID not in d:
                raise Exception('No chain %s in %s' % (template.chainID,
                                                       dssppath))
            seq, secstr, disulifids = d[template.chainID]
            return seq, secstr
        else:
            obj = self.yasara.LoadPDB(template.pdbac, download='yes')[0]
            sequence = self.yasara.SequenceMol('obj %i and mol %s'
                                               % (obj, template.chainID))
            secstr = self.yasara.SecStrRes('obj %i and mol %s'
                                           % (obj, template.chainID))
            self.yasara.DelObj(obj)
            if 'X' in secstr:
                raise Exception(
                        "yasara reported secondary structure %s for %s_%s"
                        % (secstr, template.pdbac, template.chainID))
            return seq, secstr

secstr = SecondaryStructureProvider()
