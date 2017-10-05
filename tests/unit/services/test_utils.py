from nose.tools import eq_, ok_
from mock import patch

from hommod_rest.services.modelutils import TemplateID, BlastAlignment
from hommod_rest.services.utils import extract_template_id, blast_models


@patch('hommod_rest.services.utils.extract_info')
def test_extract_template_id(mock_extract_info):
    mock_extract_info.return_value = {
        'template': '1xxx_A'
    }
    rv = extract_template_id('test/path/not/important')
    eq_(rv, TemplateID('1xxx', 'A'))


@patch('hommod_rest.services.blast.blaster.blast_models')
@patch('os.path.isfile')
@patch('hommod_rest.services.utils.extract_info')
def test_blast_models(mock_extract_info, mock_isfile, mock_blast_models):

    alignment = BlastAlignment(
        28, 193,
"APPRLICDSRVLERYLLEAKEAENITTGCAEHCSLNENITVPDTKVNFYAWKRMEVGQQA" + \
"VEVWQGLALLSEAVLRGQALLVNSSQPWEPLQLHVDKAVSGLRSLTTLLRALGAQKEAIS" + \
"PPDAASAAPLRTITADTFRKLFRVYSNFLRGKLKLYTGEACRTGDR",
        1, 166,
"APPRLICDSRVLERYLLEAKEAENITTGCAEHCSLNENITVPDTKVNFYAWKRMEVGQQA" + \
"VEVWQGLALLSEAVLRGQALLVNSSQPWEPLQLHVDKAVSGLRSLTTLLRALGAQKEAIS" + \
"PPDAASAAPLRTITADTFRKLFRVYSNFLRGKLKLYTGEACRTGDR"
    )
    mock_blast_models.return_value = {
        '23c094b6-36fa-a580-2240-fb1cf90738d3_HUMAN_22-129|A': [alignment]
    }

    mock_isfile.return_value = True

    mock_extract_info.return_value = {
        'template': '1xxx_A'
    }

    seq = \
"MGVHECPAWLWLLLSLLSLPLGLPVLGAPPRLICDSRVLERYLLEAKEAENITTGCAEHC" + \
"SLNENITVPDTKVNFYAWKRMEVGQQAVEVWQGLALLSEAVLRGQALLVNSSQPWEPLQL" + \
"HVDKAVSGLRSLTTLLRALGAQKEAISPPDAASAAPLRTITADTFRKLFRVYSNFLRGKL" + \
"KLYTGEACRTGDR"

    approved = blast_models(seq, 'HUMAN', 137, None)
    ok_(len(approved) > 0)
