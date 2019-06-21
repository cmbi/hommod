from nose.tools import with_setup, ok_

from hommod import default_settings as settings
from hommod.controllers.storage import model_storage
from hommod.controllers.clustal import clustal_aligner


def setup():
    model_storage.model_dir = "data"
    clustal_aligner.clustalw_exe = settings.CLUSTALW_EXE


def teardown():
    pass


@with_setup(setup, teardown)
def test_covered():
    sequence = "MELLTFRDVAIEFSPEEWKCLDPDQQNLYRDVMLENYRNLVSLGVAISNPDLVTCLEQRKEPYNVKIHKIVARPPAMCSHFTQDHWPVQGIEDSFHKLILRRYEKCGHDNLQLRKGCKSLNECKLQKGGYNEFNECLSTTQSKILQCKASVKVVSKFSNSNKRKTRHTGEKHFKECGKSFQKFSHLTQHKVIHAGEKPYTCEECGKAFKWSLIFNEHKRIHTGEKPFTCEECGSIFTTSSHFAKHKIIHTGEKPYKCEECGKAFNRFTTLTKHKRIHAGEKPITCEECRKIFTSSSNFAKHKRIHTGEKPYKCEECGKAFNRSTTLTKHKRIHTGEKPYTCEECGKAFRQSSKLNEHKKVHTGERPYKCDECGKAFGRSRVLNEHKKIHTGEKPYKCEECGKAFRRSTDRSQHKKIHSADKPYKCKECDKAFKQFSLLSQHKKIHTVDKPYKCKDCDKAFKRFSHLNKHKKIHT"
    position = 385

    ok_(not model_storage.model_covers("tests/unit/data/zfn.tgz", sequence, position))
