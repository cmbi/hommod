from nose.tools import ok_
from mock import patch

from hommod.controllers.clustal import clustal_aligner
from hommod.controllers.method import select_best_model
from hommod.models.align import Alignment
from hommod.models.aminoacid import AminoAcid


@patch("hommod.controllers.storage.model_storage.get_model_name_from_path")
@patch("hommod.controllers.storage.model_storage.get_sequence_id_from_name")
@patch("hommod.controllers.storage.model_storage.extract_alignments")
@patch("hommod.controllers.storage.model_storage.extract_selected_targets")
@patch("hommod.controllers.storage.model_storage.extract_model")
def test_select_best_model(mock_extract_model,
                           mock_extract_selected_targets,
                           mock_extract_alignments,
                           mock_get_sequence_id_from_name,
                           mock_get_model_name_from_path):

    from hommod import default_settings as settings
    clustal_aligner.clustalw_exe = settings.CLUSTALW_EXE

    target_id = 'target'

    paths = ['0', '1', '2']

    alignments = [
        Alignment({target_id: "ATATATATATAT", '1xxx': "ATATATVTATAT"}),
        Alignment({target_id: "ATATATATATAT", '2xxx': "AVAV--VTVTAY"}),
        Alignment({target_id: "ATATATATATAT", '3xxx': "ATATA-ATATAT"})
    ]

    pdb_str = """
SEQRES   1 A  154  ALA THR ALA THR ALA THR ALA THR ALA THR ALA THR
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C
ATOM      2  CA  ALA A   2       0.000   0.000   0.000  1.00 20.00           C
ATOM      3  CA  ALA A   3       0.000   0.000   0.000  1.00 20.00           C
"""
    mock_extract_model.return_value = pdb_str

    mock_get_model_name_from_path.return_value = 'fake-name'
    mock_get_sequence_id_from_name.return_value = target_id

    def extract_alignments(path):
        return [alignments[int(path)]]

    mock_extract_alignments.side_effect = extract_alignments

    def extract_selected_targets(path):
        return {'A': target_id}

    mock_extract_selected_targets.side_effect = extract_selected_targets

    best = select_best_model(paths, 'ATATATATATAT', 2)
    ok_(best is not None)
