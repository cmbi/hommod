import tarfile
import os
import shutil
import tempfile

from nose.tools import with_setup, ok_
from celery import Celery

from hommod.default_settings import (KMAD_EXE, BLASTP_EXE, MODEL_DIR, INTERPRO_URL,
                                     FORBIDDEN_INTERPRO_DOMAINS, TEMPLATE_BLAST_DATABANK,
                                     DOMAIN_MIN_PERCENTAGE_COVERAGE, DOMAIN_MAX_MERGE_DISTANCE,
                                     DSSP_DIR, SIMILAR_RANGES_MIN_OVERLAP_PERCENTAGE,
                                     SIMILAR_RANGES_MAX_LENGTH_DIFFERENCE_PERCENTAGE,
                                     BLACKLIST_FILE_PATH, HIGHLY_HOMOLOGOUS_PERCENTAGE_IDENTITY,
                                     CELERY_RESULT_BACKEND, CELERY_BROKER_URL, YASARA_DIR,
                                     CLUSTALW_EXE, UNIPROT_BLAST_DATABANK, UNIPROT_URL,
                                     CACHE_REDIS_HOST, CACHE_REDIS_PORT, CACHE_REDIS_DB,
                                     CACHE_EXPIRATION_TIME, CACHE_LOCK_TIMEOUT)
from hommod.controllers.kmad import kmad_aligner
from hommod.controllers.clustal import clustal_aligner
from hommod.controllers.model import modeler
from hommod.controllers.soup import soup
from hommod.services.interpro import interpro
from hommod.controllers.domain import domain_aligner
from hommod.controllers.blast import blaster
from hommod.controllers.blacklist import blacklister
from hommod.services.dssp import dssp
from hommod.services.uniprot import uniprot
from hommod.models.template import TemplateID
from hommod.controllers.storage import model_storage
from hommod.services.helpers.cache import cache_manager as cm


def setup():
    uniprot.url = UNIPROT_URL
    interpro.url = INTERPRO_URL
    kmad_aligner.kmad_exe = KMAD_EXE
    clustal_aligner.clustalw_exe = CLUSTALW_EXE
    model_storage.model_dir = tempfile.mkdtemp()
    soup.yasara_dir = YASARA_DIR
    modeler.uniprot_databank = UNIPROT_BLAST_DATABANK
    domain_aligner.forbidden_interpro_domains = FORBIDDEN_INTERPRO_DOMAINS
    domain_aligner.similar_ranges_min_overlap_percentage = SIMILAR_RANGES_MIN_OVERLAP_PERCENTAGE
    domain_aligner.similar_ranges_max_length_difference_percentage = SIMILAR_RANGES_MAX_LENGTH_DIFFERENCE_PERCENTAGE
    domain_aligner.min_percentage_coverage = DOMAIN_MIN_PERCENTAGE_COVERAGE
    domain_aligner.template_blast_databank = TEMPLATE_BLAST_DATABANK
    domain_aligner.max_merge_distance = DOMAIN_MAX_MERGE_DISTANCE
    domain_aligner.highly_homologous_percentage_identity = HIGHLY_HOMOLOGOUS_PERCENTAGE_IDENTITY
    blaster.blastp_exe = BLASTP_EXE
    dssp.dssp_dir = DSSP_DIR
    blacklister.file_path = BLACKLIST_FILE_PATH
    cm.redis_hostname = CACHE_REDIS_HOST
    cm.redis_port = CACHE_REDIS_PORT
    cm.redis_db = CACHE_REDIS_DB
    cm.expiration_time = CACHE_EXPIRATION_TIME
    cm.lock_timeout = CACHE_LOCK_TIMEOUT

    celery_app = Celery(__name__,
                        backend=CELERY_RESULT_BACKEND,
                        broker=CELERY_BROKER_URL)
    celery_app.conf.update({'TESTING': True, 'CELERY_ALWAYS_EAGER': True})


def end():
    shutil.rmtree(model_storage.model_dir)


@with_setup(setup, end)
def test_create_model_crambin():
    from hommod.tasks import create_model

    path = create_model("VTCCPSIVARSNFNVCRLPGTPQALCATYTGCIIIPGATCPGDFAN", "CRAAB")
    ok_(path is not None)

    name = os.path.splitext(os.path.basename(path))[0]
    pdb_name = os.path.join(name, 'target.pdb')
    with tarfile.open(path) as tf:
        ok_(pdb_name in tf.getnames())


@with_setup(setup, end)
def test_create_model_cox():
    from hommod.tasks import create_model

    path = create_model(
"MLATRVFSLVGKRAISTSVCVRAHESVVKSEDFSLPAYMDRRDHPLPEVAHVKHLSASQKALKEKEKASWSS" +
"LSMDEKVELYRIKFKESFAEMNRGSNEWKTVVGGAMFFIGFTALVIMWQKHYVYGPLPQSFDKEWVAKQTKR" +
"MLDMKVNPIQGLASKWDYEKNEWKK", "HUMAN", None, TemplateID('2y69', 'D'))
    ok_(path is not None)

    name = os.path.splitext(os.path.basename(path))[0]
    pdb_name = os.path.join(name, 'target.pdb')
    with tarfile.open(path) as tf:
        ok_(pdb_name in tf.getnames())

    alignments = model_storage.extract_alignments(path)
    ok_(len(alignments) >= 13)


@with_setup(setup, end)
def test_create_model_hydrolase():
    from hommod.tasks import create_model

    path = create_model(
"MSRIEKMSILGVRSFGIEDKDKQIITFFSPLTILVGPNGAGKTTIIECLKYICTGDFPPGTKGNTFVHDPKV"
"AQETDVRAQIRLQFRDVNGELIAVQRSMVCTQKSKKTEFKTLEGVITRTKHGEKVSLSSKCAEIDREMISSL"
"GVSKAVLNNVIFCHQEDSNWPLSEGKALKQKFDEIFSATRYIKALETLRQVRQTQGQKVKEYQMELKYLKQY"
"KEKACEIRDQITSKEAQLTSSKEIVKSYENELDPLKNRLKEIEHNLSKIMKLDNEIKALDSRKKQMEKDNSE"
"LEEKMEKVFQGTDEQLNDLYHNHQRTVREKERKLVDCHRELEKLNKESRLLNQEKSELLVEQGRLQLQADRH"
"QEHIRARDSLIQSLATQLELDGFERGPFSERQIKNFHKLVRERQEGEAKTANQLMNDFAEKETLKQKQIDEI"
"RDKKTGLGRIIELKSEILSKKQNELKNVKYELQQLEGSSDRILELDQELIKAERELSKAEKNSNVETLKMEV"
"ISLQNEKADLDRTLRKLDQEMEQLNHHTTTRTQMEMLTKDKADKDEQIRKIKSRHSDELTSLLGYFPNKKQL"
"EDWLHSKSKEINQTRDRLAKLNKELASSEQNKNHINNELKRKEEQLSSYEDKLFDVCGSQDFESDLDRLKEE"
"IEKSSKQRAMLAGATAVYSQFITQLTDENQSCCPVCQRVFQTEAELQEVISDLQSKLRLAPDKLKSTESELK"
"KKEKRRDEMLGLVPMRQSIIDLKEKEIPELRNKLQNVNRDIQRLKNDIEEQETLLGTIMPEEESAKVCLTDV"
"TIMERFQMELKDVERKIAQQAAKLQGIDLDRTVQQVNQEKQEKQHKLDTVSSKIELNRKLIQDQQEQIQHLK"
"STTNELKSEKLQISTNLQRRQQLEEQTVELSTEVQSLYREIKDAKEQVSPLETTLEKFQQEKEELINKKNTS"
"NKIAQDKLNDIKEKVKNIHGYMKDIENYIQDGKDDYKKQKETELNKVIAQLSECEKHKEKINEDMRLMRQDI"
"DTQKIQERWLQDNLTLRKRNEELKEVEEERKQHLKEMGQMQVLQMKSEHQKLEENIDNIKRNHNLALGRQKG"
"YEEEIIHFKKELREPQFRDAEEKYREMMIVMRTTELVNKDLDIYYKTLDQAIMKFHSMKMEEINKIIRDLWR"
"STYRGQDIEYIEIRSDADENVSASDKRRNYNYRVVMLKGDTALDMRGRCSAGQKVLASLIIRLALAETFCLN"
"CGIIALDEPTTNLDRENIESLAHALVEIIKSRSQQRNFQLLVITHDEDFVELLGRSEYVEKFYRIKKNIDQC"
"SEIVKCSVSSLGFNVH", "HUMAN", 1184, None)
    ok_(path is not None)

    ok_('5GOX' not in path.upper())


@with_setup(setup, end)
def test_create_model_5MHF():
    from hommod.tasks import create_model

    path = create_model(
"MARGERRRRAVPAEGVRTAERAARGGPGRRDGRGGGPRSTAGGVALAVVVLSLALGMSGRWVLAWYRARRAV"
"TLHSAPPVLPADSSSPAVAPDLFWGTYRPHVYFGMKTRSPKPLLTGLMWAQQGTTPGTPKLRHTCEQGDGVG"
"PYGWEFHDGLSFGRQHIQDGALRLTTEFVKRPGGQHGGDWSWRVTVEPQDSGTSALPLVSLFFYVVTDGKEV"
"LLPEVGAKGQLKFISGHTSELGDFRFTLLPPTSPGDTAPKYGSYNVFWTSNPGLPLLTEMVKSRLNSWFQHR"
"PPGAPPERYLGLPGSLKWEDRGPSGQGQGQFLIQQVTLKIPISIEFVFESGSAQAGGNQALPRLAGSLLTQA"
"LESHAEGFRERFEKTFQLKEKGLSSGEQVLGQAALSGLLGGIGYFYGQGLVLPDIGVEGSEQKVDPALFPPV"
"PLFTAVPSRSFFPRGFLWDEGFHQLVVQRWDPSLTREALGHWLGLLNADGWIGREQILGDEARARVPPEFLV"
"QRAVHANPPTLLLPVAHMLEVGDPDDLAFLRKALPRLHAWFSWLHQSQAGPLPLSYRWRGRDPALPTLLNPK"
"TLPSGLDDYPRASHPSVTERHLDLRCWVALGARVLTRLAEHLGEAEVAAELGPLAASLEAAESLDELHWAPE"
"LGVFADFGNHTKAVQLKPRPPQGLVRVVGRPQPQLQYVDALGYVSLFPLLLRLLDPTSSRLGPLLDILADSR"
"HLWSPFGLRSLAASSSFYGQRNSEHDPPYWRGAVWLNVNYLALGALHHYGHLEGPHQARAAKLHGELRANVV"
"GNVWRQYQATGFLWEQYSDRDGRGMGCRPFHGWTSLVLLAMAEDY", "HUMAN", 100,
                        TemplateID('5MHF', 'D'))
    ok_(path is not None)

    contents = model_storage.extract_model(path)
    ok_(contents is not None)
