from nose.tools import with_setup, ok_

from hommod.services.uniprot import uniprot


def setup():
    import hommod.default_settings as settings
    uniprot.url = settings.UNIPROT_URL


def teardown():
    pass


@with_setup(setup, teardown)
def test_get_sequence():
    sequence = uniprot.get_sequence('P01542')
    ok_(len(sequence) > 0)
