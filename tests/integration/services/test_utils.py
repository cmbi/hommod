from nose.tools import ok_, with_setup

from hommod_rest.services.utils import blast_models
from hommod_rest.services.blast import blaster
import hommod_rest.default_settings as settings


def setup():
    blaster.blastp_exe = settings.BLASTP
    blaster.templates_db = settings.TEMPLATESDB
    blaster.uniprot_db = settings.UNIPROTDB
    blaster.models_db = settings.MODELSDB

def end():
    pass

@with_setup(setup, end)
def test_blast_models():
    seq = \
"MGVHECPAWLWLLLSLLSLPLGLPVLGAPPRLICDSRVLERYLLEAKEAENITTGCAEHC" + \
"SLNENITVPDTKVNFYAWKRMEVGQQAVEVWQGLALLSEAVLRGQALLVNSSQPWEPLQL" + \
"HVDKAVSGLRSLTTLLRALGAQKEAISPPDAASAAPLRTITADTFRKLFRVYSNFLRGKL" + \
"KLYTGEACRTGDR"

    approved = blast_models(seq, 'HUMAN', 137, None)
    ok_(len(approved) > 0)
