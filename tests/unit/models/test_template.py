from nose.tools import eq_, ok_

from hommod.models.template import TemplateID


def test_compare_template_id():
    template_id_1 = TemplateID('2ypd', 'A')
    template_id_2 = TemplateID('2ypd', 'A')
    template_id_3 = TemplateID('2ypd', 'B')

    eq_(template_id_1, template_id_2)
    ok_(template_id_1 != template_id_3)
    ok_(not (template_id_1 != template_id_2))
