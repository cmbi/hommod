from nose.tools import eq_
from mock import patch

from hommod_rest.services.modelutils import TemplateID
from hommod_rest.services.utils import extract_template_id


@patch('hommod_rest.services.utils.extract_info')
def test_extract_template_id(mock_extract_info):
    mock_extract_info.return_value = {
        'template': '1xxx_A'
    }
    rv = extract_template_id('test/path/not/important')
    eq_(rv, TemplateID('1xxx', 'A'))
