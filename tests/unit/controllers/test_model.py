import tempfile
import os

from nose.tools import eq_
from mock import patch

from hommod.controllers.model import modeler
from hommod.controllers.storage import model_storage
from hommod.models.align import TargetTemplateAlignment


@patch("tarfile.open")
def test_selected_targets_format(mock_open):

    targets = {'a': 'coffee',
               'b': 'tea'}

    alignments = {}
    for chain_id in targets:
        alignments[chain_id] = TargetTemplateAlignment('', '')
        alignments[chain_id].target_id = targets[chain_id]

    path = tempfile.mktemp()

    class FakeTarFile:
        def __enter__(self, *args, **kwargs):
            return self

        def __exit__(self, *args, **kwargs):
            pass

        def extractfile(self, *args, **kwargs):
            return open(path, 'rb')

    mock_open.return_value = FakeTarFile()

    try:
        modeler._write_selected_targets(alignments, path)
        parsed = model_storage.extract_selected_targets('no.tar.gz')
    finally:
        if os.path.isfile(path):
            os.remove(path)

    eq_(parsed, targets)
