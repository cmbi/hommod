import tarfile
import os
from glob import glob
import shutil
import tempfile
import logging

from mock import patch
from nose.tools import with_setup, ok_, raises
from celery import Celery

from hommod.default_settings import (KMAD_EXE, BLASTP_EXE, MODEL_DIR, INTERPRO_URL,
                                     FORBIDDEN_INTERPRO_DOMAINS, TEMPLATE_BLAST_DATABANK,
                                     DOMAIN_MIN_PERCENTAGE_COVERAGE,
                                     DSSP_DIR, SIMILAR_RANGES_MIN_OVERLAP_PERCENTAGE,
                                     SIMILAR_RANGES_MAX_LENGTH_DIFFERENCE_PERCENTAGE,
                                     BLACKLIST_FILE_PATH, HIGHLY_HOMOLOGOUS_PERCENTAGE_IDENTITY,
                                     YASARA_DIR,
                                     CLUSTALW_EXE, UNIPROT_BLAST_DATABANK, SPROT_FASTA, TREMBL_FASTA,
                                     CACHE_REDIS_HOST, CACHE_REDIS_PORT, CACHE_REDIS_DB,
                                     CACHE_EXPIRATION_TIME, CACHE_LOCK_TIMEOUT)
from hommod.controllers.kmad import kmad_aligner
from hommod.controllers.clustal import clustal_aligner
from hommod.controllers.model import modeler
from hommod.services.interpro import interpro
from hommod.controllers.domain import domain_aligner
from hommod.controllers.blast import blaster
from hommod.controllers.blacklist import blacklister
from hommod.services.dssp import dssp
from hommod.controllers.uniprot import uniprot
from hommod.models.template import TemplateID
from hommod.controllers.storage import model_storage
from hommod.controllers.method import select_best_model
from hommod.services.helpers.cache import cache_manager as cm
from hommod.models.error import ModelRunError


_log = logging.getLogger(__name__)


def setup():
    uniprot.fasta_paths = [SPROT_FASTA, TREMBL_FASTA]
    interpro.url = INTERPRO_URL
    interpro.email = "test.integration@radboudumc.nl"
    kmad_aligner.kmad_exe = KMAD_EXE
    clustal_aligner.clustalw_exe = CLUSTALW_EXE
    model_storage.model_dir = tempfile.mkdtemp()
    modeler.yasara_dir = YASARA_DIR
    modeler.uniprot_databank = UNIPROT_BLAST_DATABANK
    domain_aligner.forbidden_interpro_domains = FORBIDDEN_INTERPRO_DOMAINS
    domain_aligner.similar_ranges_min_overlap_percentage = SIMILAR_RANGES_MIN_OVERLAP_PERCENTAGE
    domain_aligner.similar_ranges_max_length_difference_percentage = SIMILAR_RANGES_MAX_LENGTH_DIFFERENCE_PERCENTAGE
    domain_aligner.min_percentage_coverage = DOMAIN_MIN_PERCENTAGE_COVERAGE
    domain_aligner.template_blast_databank = TEMPLATE_BLAST_DATABANK
    domain_aligner.highly_homologous_percentage_identity = HIGHLY_HOMOLOGOUS_PERCENTAGE_IDENTITY
    blaster.blastp_exe = BLASTP_EXE
    dssp.dssp_dir = DSSP_DIR
    blacklister.file_path = BLACKLIST_FILE_PATH
    cm.disable()

    celery_app = Celery(__name__)
    celery_app.conf.update({'TESTING': True, 'CELERY_ALWAYS_EAGER': True})


def end():
    cm.enable()
    shutil.rmtree(model_storage.model_dir)


@with_setup(setup, end)
def test_create_model_1pfx():
    from hommod.tasks import create_model

    path = create_model("MQRVNMIMAESPGLITICLLGYLLSAECTVFLDHENANKILNRPKRYNSGKLEEFVQGNLERECMEEKCSFEEAREVFENTERTTEFWKQYVDGDQCESNPCLNGGSCKDDINSYECWCPFGFEGKNCELDVTCNIKNGRCEQFCKNSADNKVVCSCTEGYRLAENQKSCEPAVPFPCGRVSVSQTSKLTRAETVFPDVDYVNSTEAETILDNITQSTQSFNDFTRVVGGEDAKPGQFPWQVVLNGKVDAFCGGSIVNEKWIVTAAHCVETGVKITVVAGEHNIEETEHTEQKRNVIRIIPHHNYNAAINKYNHDIALLELDEPLVLNSYVTPICIADKEYTNIFLKFGSGYVSGWGRVFHKGRSALVLQYLRVPLVDRATCLRSTKFTIYNNMFCAGFHEGGRDSCQGDSGGPHVTEVEGTSFLTGIISWGEECAMKGKYGIYTKVSRYVNWIKEKTKLT",
                        "HUMAN", 180,  TemplateID('1pfx', 'L'))
    ok_(path is not None)


@with_setup(setup, end)
def test_create_model_neurexin():
    from hommod.tasks import create_model

    path = create_model("MWLLALCLVGLAGAQRGGGGPGGGAPGGPGLGLGSLGEERFPVVNTAYGRVRGVRRELNNEILGPVVQFLGVPYATPPLGARRFQPPEAPASWPGVRNATTLPPACPQNLHGALPAIMLPVWFTDNLEAAATYVQNQSEDCLYLNLYVPTEDGPLTKKRDEATLNPPDTDIRDPGKKPVMLFLHGGSYMEGTGNMFDGSVLAAYGNVIVATLNYRLGVLGFLSTGDQAAKGNYGLLDQIQALRWLSENIAHFGGDPERITIFGSGAGASCVNLLILSHHSEGLFQKAIAQSGTAISSWSVNYQPLKYTRLLAAKVGCDREDSAEAVECLRRKPSRELVDQDVQPARYHIAFGPVVDGDVVPDDPEILMQQGEFLNYDMLIGVNQGEGLKFVEDSAESEDGVSASAFDFTVSNFVDNLYGYPEGKDVLRETIKFMYTDWADRDNGEMRRKTLLALFTDHQWVAPAVATAKLHADYQSPVYFYTFYHHCQAEGRPEWADAAHGDELPYVFGVPMVGATDLFPCNFSKNDVMLSAVVMTYWTNFAKTGDPNQPVPQDTKFIHTKPNRFEEVVWSKFNSKEKQYLHIGLKPRVRDNYRANKVAFWLELVPHLHNLHTELFTTTTRLPPYATRWPPRPPAGAPGTRRPPPPATLPPEPEPEPGPRAYDRFPGDSRDYSTELSVTVAVGASLLFLNILAFAALYYKRDRRQELRCRRLSPPGGSGSGVPGGGPLLPAAGRELPPEEELVSLQLKRGGGVGADPAEALRPACPPDYTLALRRAPDDVPLLAPGALTLLPSGLGPPPPPPPPSLHPFGPFPPPPPTATSHNNTLPHPHSTTRV",
                        "HUMAN", 81, TemplateID('3vkf', 'A'))
    ok_(path is not None)


@with_setup(setup, end)
def test_create_model_crambin():
    from hommod.tasks import create_model

    path = create_model("VTCCPSIVARSNFNVCRLPGTPQALCATYTGCIIIPGATCPGDFAN", "CRAAB")
    ok_(path is not None)

    name = os.path.splitext(os.path.basename(path))[0]
    pdb_name = os.path.join(name, 'target.pdb')
    log_name = os.path.join(name, 'model.log')
    with tarfile.open(path) as tf:
        ok_(pdb_name in tf.getnames())
        ok_(log_name in tf.getnames())


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


@with_setup(setup, end)
def test_create_model_5GJV():
    from hommod.tasks import create_model

    path = create_model(
"MEPSSPQDEGLRKKQPKKPVPEILPRPPRALFCLTLENPLRKACISIVEWKPFETIILLTIFANCVALAVYL" +
"PMPEDDNNSLNLGLEKLEYFFLIVFSIEAAMKIIAYGFLFHQDAYLRSGWNVLDFTIVFLGVFTVILEQVNV" +
"IQSHTAPMSSKGAGLDVKALRAFRVLRPLRLVSGVPSLQVVLNSIFKAMLPLFHIALLVLFMVIIYAIIGLE" +
"LFKGKMHKTCYFIGTDIVATVENEEPSPCARTGSGRRCTINGSECRGGWPGPNHGITHFDNFGFSMLTVYQC" +
"ITMEGWTDVLYWVNDAIGNEWPWIYFVTLILLGSFFILNLVLGVLSGEFTKEREKAKSRGTFQKLREKQQLD" +
"EDLRGYMSWITQGEVMDVEDFREGKLSLDEGGSDTESLYEIAGLNKIIQFIRHWRQWNRIFRWKCHDIVKSK" +
"VFYWLVILIVALNTLSIASEHHNQPLWLTRLQDIANRVLLSLFTTEMLMKMYGLGLRQYFMSIFNRFDCFVV" +
"CSGILEILLVESGAMTPLGISVLRCIRLLRIFKITKYWTSLSNLVASLLNSIRSIASLLLLLFLFIVIFALL" +
"GMQLFGGRYDFEDTEVRRSNFDNFPQALISVFQVLTGEDWTSMMYNGIMAYGGPSYPGMLVCIYFIILFVCG" +
"NYILLNVFLAIAVDNLAEAESLTSAQKAKAEEKKRRKMSKGLPDKSEEEKSTMAKKLEQKPKGEGIPTTAKL" +
"KIDEFESNVNEVKDPYPSADFPGDDEEDEPEIPLSPRPRPLAELQLKEKAVPIPEASSFFIFSPTNKIRVLC" +
"HRIVNATWFTNFILLFILLSSAALAAEDPIRADSMRNQILKHFDIGFTSVFTVEIVLKMTTYGAFLHKGSFC" +
"RNYFNMLDLLVVAVSLISMGLESSAISVVKILRVLRVLRPLRAINRAKGLKHVVQCMFVAISTIGNIVLVTT" +
"LLQFMFACIGVQLFKGKFFRCTDLSKMTEEECRGYYYVYKDGDPMQIELRHREWVHSDFHFDNVLSAMMSLF" +
"TVSTFEGWPQLLYKAIDSNAEDVGPIYNNRVEMAIFFIIYIILIAFFMMNIFVGFVIVTFQEQGETEYKNCE" +
"LDKNQRQCVQYALKARPLRCYIPKNPYQYQVWYIVTSSYFEYLMFALIMLNTICLGMQHYNQSEQMNHISDI" +
"LNVAFTIIFTLEMILKLMAFKARGYFGDPWNVFDFLIVIGSIIDVILSEIDTFLASSGGLYCLGGGCGNVDP" +
"DESARISSAFFRLFRVMRLIKLLSRAEGVRTLLWTFIKSFQALPYVALLIVMLFFIYAVIGMQMFGKIALVD" +
"GTQINRNNNFQTFPQAVLLLFRCATGEAWQEILLACSYGKLCDPESDYAPGEEYTCGTNFAYYYFISFYMLC" +
"AFLVINLFVAVIMDNFDYLTRDWSILGPHHLDEFKAIWAEYDPEAKGRIKHLDVVTLLRRIQPPLGFGKFCP" +
"HRVACKRLVGMNMPLNSDGTVTFNATLFALVRTALKIKTEGNFEQANEELRAIIKKIWKRTSMKLLDQVIPP" +
"IGDDEVTVGKFYATFLIQEHFRKFMKRQEEYYGYRPKKDIVQIQAGLRTIEEEAAPEICRTVSGDLAAEEEL" +
"ERAMVEAAMEEGIFRRTGGLFGQVDNFLERTNSLPPVMANQRPLQFAEIEMEEMESPVFLEDFPQDPRTNPL" +
"ARANTNNANANVAYGNSNHSNSHVFSSVHYEREFPEETETPATRGRALGQPCRVLGPHSKPCVEMLKGLLTQ" +
"RAMPRGQAPPAPCQCPRVESSMPEDRKSSTPGSLHEETPHSRSTRENTSRCSAPATALLIQKALVRGGLGTL" +
"AADANFIMATGQALADACQMEPEEVEIMATELLKGREAPEGMASSLGCLNLGSSLGSLDQHQGSQETLIPPR",
    "HUMAN", 434, TemplateID('5GJV', 'A'))

    ok_(path is not None)


@with_setup(setup, end)
def test_create_model_delmol():
    from hommod.tasks import create_model

    path = create_model(
"MESADFYEAEPRPPMSSHLQSPPHAPSSAAFGFPRGAGPAQPPAPPAAPEPLGGICEHETSIDISAYIDPAA"
"FNDEFLADLFQHSRQQEKAKAAVGPTGGGGGGDFDYPGAPAGPGGAVMPGGAHGPPPGYGCAAAGYLDGRLE"
"PLYERVGAPALRPLVIKQEPREEDEAKQLALAGLFPYQPPPPPPPSHPHPHPPPAHLAAPHLQFQIAHCGQT"
"TMHLQPGHPTPPPTPVPSPHPAPALGAAGLPGPGSALKGLGAAHPDLRASGGSGAGKAKKSVDKNSNEYRVR"
"RERNNIAVRKSRDKAKQRNVETQQKVLELTSDNDRLRKRVEQLSRELDTLRGIFRQLPESSLVKAMGNCA",
    "HUMAN", 317, TemplateID('1H88', 'B'))

    ok_(path is not None)


@with_setup(setup, end)
def test_create_model_cadherin():
    from hommod.tasks import create_model

    path = create_model(
           "MTIHQFLLLFLFWVCLPHFCSPEIMFRRTPVPQQRILSSRVPRSDGKILHRQKRGWMWNQ" +
           "FFLLEEYTGSDYQYVGKLHSDQDKGDGSLKYILSGDGAGTLFIIDEKTGDIHATRRIDRE" +
           "EKAFYTLRAQAINRRTLRPVEPESEFVIKIHDINDNEPTFPEEIYTASVPEMSVVGTSVV" +
           "QVTATDADDPSYGNSARVIYSILQGQPYFSVEPETGIIRTALPNMNRENREQYQVVIQAK" +
           "DMGGQMGGLSGTTTVNITLTDVNDNPPRFPQNTIHLRVLESSPVGTAIGSVKATDADTGK" +
           "NAEVEYRIIDGDGTDMFDIVTEKDTQEGIITVKKPLDYESRRLYTLKVEAENTHVDPRFY" +
           "YLGPFKDTTIVKISIEDVDEPPVFSRSSYLFEVHEDIEVGTIIGTVMARDPDSISSPIRF" +
           "SLDRHTDLDRIFNIHSGNGSLYTSKPLDRELSQWHNLTVIAAEINNPKETTRVAVFVRIL" +
           "DVNDNAPQFAVFYDTFVCENARPGQLIQTISAVDKDDPLGGQKFFFSLAAVNPNFTVQDN" +
           "EDNTARILTRKNGFNRHEISTYLLPVVISDNDYPIQSSTGTLTIRVCACDSQGNMQSCSA" +
           "EALLLPAGLSTGALIAILLCIIILLVIVVLFAALKRQRKKEPLILSKEDIRDNIVSYNDE" +
           "GGGEEDTQAFDIGTLRNPAAIEEKKLRRDIIPETLFIPRRTPTAPDNTDVRDFINERLKE" +
           "HDLDPTAPPYDSLATYAYEGNDSIAESLSSLESGTTEGDQNYDYLREWGPRFNKLAEMYG" +
           "GGESDKDS", "HUMAN", 540)

    ok_(path is not None)


@with_setup(setup, end)
def test_create_model_LSHR_extracellular():
    from hommod.tasks import create_model

    path = create_model(
"MKQRFSALQLLKLLLLLQPPLPRALREALCPEPCNCVPDGALRCPGPTAGLTRLSLAYLPVKVIPSQAFRGLNEVIKIEISQIDSLERIEANAFDNLLNLSEILIQNTKNLRYIEPGAFINLPRLKYLSICNTGIRKFPDVTKVFSSESNFILEICDNLHITTIPGNAFQGMNNESVTLKLYGNGFEEVQSHAFNGTTLTSLELKENVHLEKMHNGAFRGATGPKTLDISSTKLQALPSYGLESIQRLIATSSYSLKKLPSRETFVNLLEATLTYPSHCCAFRNLPTKEQNFSHSISENFSKQCESTVRKVNNKTLYSSMLAESELSGWDYEYGFCLPKTPRCAPEPDAFNPCEDIMGYDFLRVLIWLINILAIMGNMTVLFVLLTSRYKLTVPRFLMCNLSFADFCMGLYLLLIASVDSQTKGQYYNHAIDWQTGSGCSTAGFFTVFASELSVYTLTVITLERWHTITYAIHLDQKLRLRHAILIMLGGWLFSSLIAMLPLVGVSNYMKVSICFPMDVETTLSQVYILTILILNVVAFFIICACYIKIYFAVRNPELMATNKDTKIAKKMAILIFTDFTCMAPISFFAISAAFKVPLITVTNSKVLLVLFYPINSCANPFLYAIFTKTFQRDFFLLLSKFGCCKRRAELYRRKDFSAYTSNCKNGFTGSNKPSQSTLKLSTLHCQGTALLDKTRYTEC",
'HUMAN', 428, TemplateID('2J4Y', 'A'))

    ok_(path is None)


@with_setup(setup, end)
def test_create_no_model_bestrophin():
    from hommod.tasks import create_model
    path = create_model(
"MTITYTSQVANARLGSFSRLLLCWRGSIYKLLYGEFLIFLLCYYIIRFIYRLALTEEQQLMFEKLTLYCDSY" +
"IQLIPISFVLGFYVTLVVTRWWNQYENLPWPDRLMSLVSGFVEGKDEQGRLLRRTLIRYANLGNVLILRSVS" +
"TAVYKRFPSAQHLVQAGFMTPAEHKQLEKLSLPHNMFWVPWVWFANLSMKAWLGGRIRDPILLQSLLNEMNT" +
"LRTQCGHLYAYDWISIPLVYTQVVTVAVYSFFLTCLVGRQFLNPAKAYPGHELDLVVPVFTFLQFFFYVGWL" +
"KVAEQLINPFGEDDDDFETNWIVDRNLQVSLLAVDEMHQDLPRMEPDMYWNKPEPQPPYTAASAQFRRASFM" +
"GSTFNISLNKEEMEFQPNQEDEEDAHAGIIGRFLGLQSHDHHPPRANSRTKLLWPKRESLLHEGLPKNHKAA" +
"KQNVRGQEDNKAWKLKAVDAFKSAPLYQRPGYYSAPQTPLSPTPMFFPLEPSAPSKLHSVTGIDTKDKSLKT" +
"VSSGAKKSFELLSESDGALMEHPEVSQVRRKTVEFNLTDMPEIPENHLKEPLEQSPTNIHTTLKDHMDPYWA" +
"LENRDEAHS", "HUMAN", 325, TemplateID('4RDQ', 'A'))

    ok_(path is None)


@with_setup(setup, end)
def test_create_model_gprotein():
    from hommod.tasks import create_model
    path = create_model(
"MTLESIMACCLSEEAKEARRINDEIERQLRRDKRDARRELKLLLLGTGESGKSTFIKQMRIIHGSGYSDEDKRGFTKLVYQNIFTAMQAMIRAMDTLKIPYKYEHNKAHAQLVREVDVEKVSAFENPYVDAIKSLWNDPGIQECYDRRREYQLSDSTKYYLNDLDRVADPAYLPTQQDVLRVRVPTTGIIEYPFDLQSVIFRMVDVGGQRSERRKWIHCFENVTSIMFLVALSEYDQVLVESDNENRMEESKALFRTIITYPWFQNSSVILFLNKKDLLEEKIMYSHLVDYFPEYDGPQRDAQAAREFILKMFVDLNPDSDKIIYSHFTCATDTENIRFVFAAVKDTILQLNLKEYNLV",
'HUMAN', 209, TemplateID('3AH8', 'A'))

    ok_(path is None)


@with_setup(setup, end)
def test_create_model_kinase():
    sequence = "METVISSDSSPAVENEHPQETPESNNSVYTSFMKSHRCYDLIPTSSKLVVFDTSLQVKKAFFALVTNGVRAAPLWDSKKQSFVVLRALSCPLGMLTITDFINILHRYYKSALVQIYELEEHKIETWREVYLQDSFKPLVCISPNASLFDAVSSLIRNKIHRLPVIDPESGNTLYILTHKRILKF"
    position = 65

    from hommod.tasks import create_model
    path = create_model(sequence, 'HUMAN', position, TemplateID('2V8Q', 'E'))

    ok_(path is not None)

    ok_(select_best_model([path], sequence, position) is not None)


@with_setup(setup, end)
def test_create_model_6F1T():
    sequence = "MANGTADVRKLFIFTTTQNYFGLMSELWDQPLLCNCLEINNFLDDGNQMLLRVQRSDAGISFSNTIEFGDTKDKVLVFFKLRPEVITDENLHDNILVSSMLESPISSLYQAVRQVFAPMLLKDQEWSRNFDPKLQNLLSELEAGLGIVLRRSDTNLTKLKFKEDDTRGILTPSDEFQFWIEQAHRGNKQISKERANYFKELFETIAREFYNLDSLSLLEVVDLVETTQDVVDDVWRQTEHDHYPESRMLHLLDIIGGSFGRFVQKKLGTLNLWEDPYYLVKESLKAGISICEQWVIVCNHLTGQVWQRYVPHPWKNEKYFPETLDKLGKRLEEVLAIRTIHEKFLYFLPASEEKIICLTRVFEPFTGLNPVQYNPYTEPLWKAAVSQYEKIIAPAEQKIAGKLKNYISEIQDSPQQLLQAFLKYKELVKRPTISKELMLERETLLARLVDSIKDFRLDFENRCRGIPGDASGPLSGKNLSEVVNSIVWVRQLELKVDDTIKIAEALLSDLPGFRCFHQSAKDLLDQLKLYEQEQFDDWSRDIQSGLSDSRSGLCIEASSRIMELDSNDGLLKVHYSDRLVILLREVRQLSALGFVIPAKIQQVANIAQKFCKQAIILKQVAHFYNSIDQQMIQSQRPMMLQSALAFEQIIKNSKAGSGGKSQITWDNPKELEGYIQKLQNAAERLATENRKLRKWHTTFCEKVVVLMNIDLLRQQQRWKDGLQELRTGLATVEAQGFQASDMHAWKQHWNHQLYKALEHQYQMGLEALNENLPEINIDLTYKQGRLQFRPPFEEIRAKYYREMKRFIGIPNQFKGVGEAGDESIFSIMIDRNASGFLTIFSKAEDLFRRLSAVLHQHKEWIVIGQVDMEALVEKHLFTVHDWEKNFKALKIKGKEVERLPSAVKVDCLNINCNPVKTVIDDLIQKLFDLLVLSLKKSIQAHLHEIDTFVTEAMEVLTIMPQSVEEIGDANLQYSKLQERKPEILPLFQEAEDKNRLLRTVAGGGLETISNLKAKWDKFELMMESHQLMIKDQIEVMKGNVKSRLQIYYQELEKFKARWDQLKPGDDVIETGQHNTLDKSAKLIKEKKIEFDDLEVTRKKLVDDCHHFRLEEPNFSLASSISKDIESCAQIWAFYEEFQQGFQEMANEDWITFRTKTYLFEEFLMNWHDRLRKVEEHSVMTVKLQSEVDKYKIVIPILKYVRGEHLSPDHWLDLFRLLGLPRGTSLEKLLFGDLLRVADTIVAKAADLKDLNSRAQGEVTIREALRELDLWGVGAVFTLIDYEDSQSRTMKLIKDWKDIVNQVGDNRCLLQSLKDSPYYKGFEDKVSIWERKLAELDEYLQNLNHIQRKWVYLEPIFGRGALPKEQTRFNRVDEDFRSIMTDIKKDNRVTTLTTHAGIRNSLLTILDQLQRCQKSLNEFLEEKRSAFPRFYFIGDDDLLEILGQSTNPSVIQSHLKKLFAGINSVCFDEKSKHITAMKSLEGEVVPFKNKVPLSNNVETWLNDLALEMKKTLEQLLKECVTTGRSSQGAVDPSLFPSQILCLAEQIKFTEDVENAIKDHSLHQIETQLVNKLEQYTNIDTSSEDPGNTESGILELKLKALILDIIHNIDVVKQLNQIQVHTTEDWAWKKQLRFYMKSDHTCCVQMVDSEFQYTYEYQGNASKLVYTPLTDKCYLTLTQAMKMGLGGNPYGPAGTGKTESVKALGGLLGRQVLVFNCDEGIDVKSMGRIFVGLVKCGAWGCFDEFNRLEESVLSAVSMQIQTIQDALKNHRTVCELLGKEVEVNSNSGIFITMNPAGKGYGGRQKLPDNLKQLFRPVAMSHPDNELIAEVILYSEGFKDAKVLSRKLVAIFNLSRELLTPQQHYDWGLRALKTVLRGSGNLLRQLNKSGTTQNANESHIVVQALRLNTMSKFTFTDCTRFDALIKDVFPGIELKEVEYDELSAALKQVFEEANYEIIPNQIKKALELYEQLCQRMGVVIVGPSGAGKSTLWRMLRAALCKTGKVVKQYTMNPKAMPRYQLLGHIDMDTREWSDGVLTNSARQVVREPQDVSSWIICDGDIDPEWIESLNSVLDDNRLLTMPSGERIQFGPNVNFVFETHDLSCASPATISRMGMIFLSDEETDLNSLIKSWLRNQPAEYRNNLENWIGDYFEKALQWVLKQNDYVVETSLVGTVMNGLSHLHGCRDHDEFIINLIRGLGGNLNMKSRLEFTKEVFHWARESPPDFHKPMDTYYDSTRGRLATYVLKKPEDLTADDFSNGLTLPVIQTPDMQRGLDYFKPWLSSDTKQPFILVGPEGCGKGMLLRYAFSQLRSTQIATVHCSAQTTSRHLLQKLSQTCMVISTNTGRVYRPKDCERLVLYLKDINLPKLDKWGTSTLVAFLQQVLTYQGFYDENLEWVGLENIQIVASMSAGGRLGRHKLTTRFTSIVRLCSIDYPEREQLQTIYGAYLEPVLHKNLKNHSIWGSSSKIYLLAGSMVQVYEQVRAKFTVDDYSHYFFTPCILTQWVLGLFRYDLEGGSSNHPLDYVLEIVAYEARRLFRDKIVGAKELHLFDIILTSVFQGDWGSDILDNMSDSFYVTWGARHNSGARAAPGQPLPPHGKPLGKLNSTDLKDVIKKGLIHYGRDNQNLDILLFHEVLEYMSRIDRVLSFPGGSLLLAGRSGVGRRTITSLVSHMHGAVLFSPKISRGYELKQFKNDLKHVLQLAGIEAQQVVLLLEDYQFVHPTFLEMINSLLSSGEVPGLYTLEELEPLLLPLKDQASQDGFFGPVFNYFTYRIQQNLHIVLIMDSANSNFMINCESNPALHKKCQVLWMEGWSNSSMKKIPEMLFSETGGGEKYNDKKRKEEKKKNSVDPDFLKSFLLIHESCKAYGATPSRYMTFLHVYSAISSSKKKELLKRQSHLQAGVSKLNEAKALVDELNRKAGEQSVLLKTKQDEADAALQMITVSMQDASEQKTELERLKHRIAEEVVKIEERKNKIDDELKEVQPLVNEAKLAVGNIKPESLSEIRSLRMPPDVIRDILEGVLRLMGIFDTSWVSMKSFLAKRGVREDIATFDARNISKEIRESVEELLFKNKGSFDPKNAKRASTAAAPLAAWVKANIQYSHVLERIHPLETEQAGLESNLKKTEDRKRKLEELLNSVGQKVSELKEKFQSRTSEAAKLEAEVSKAQETIKAAEVLINQLDREHKRWNAQVVEITEELATLPKRAQLAAAFITYLSAAPESLRKTCLEEWTKSAGLEKFDLRRFLCTESEQLIWKSEGLPSDDLSIENALVILQSRVCPFLIDPSSQATEWLKTHLKDSRLEVINQQDSNFITALELAVRFGKTLIIQEMDGVEPVLYPLLRRDLVAQGPRYVVQIGDKIIDYNEEFRLFLSTRNPNPFIPPDAASIVTEVNFTTTRSGLRGQLLALTIQHEKPDLEEQKTKLLQQEEDKKIQLAKLEESLLETLATSQGNILENKDLIESLNQTKASSALIQESLKESYKLQISLDQERDAYLPLAESASKMYFIISDLSKINNMYRFSLAAFLRLFQRALQNKQDSENTEQRIQSLISSLQHMVYEYICRCLFKADQLMFALHFVRGMHPELFQENEWDTFTGVVVGDMLRKADSQQKIRDQLPSWIDQERSWAVATLKIALPSLYQTLCFEDAALWRTYYNNSMCEQEFPSILAKKVSLFQQILVVQALRPDRLQSAMALFACKTLGLKEVSPLPLNLKRLYKETLEIEPILIIISPGADPSQELQELANAERSGECYHQVAMGQGQADLAIQMLKECARNGDWLCLKNLHLVVSWLPVLEKELNTLQPKDTFRLWLTAEVHPNFTPILLQSSLKITYESPPGLKKNLMRTYESWTPEQISKKDNTHRAHALFSLAWFHAACQERRNYIPQGWTKFYEFSLSDLRAGYNIIDRLFDGAKDVQWEFVHGLLENAIYGGRIDNYFDLRVLQSYLKQFFNSSVIDVFNQRNKKSIFPYSVSLPQSCSILDYRAVIEKIPEDDKPSFFGLPANIARSSQRMISSQVISQLRILGRSITAGSKFDREIWSNELSPVLNLWKKLNQNSNLIHQKVPPPNDRQGSPILSFIILEQFNAIRLVQSVHQSLAALSKVIRGTTLLSSEVQKLASALLNQKCPLAWQSKWEGPEDPLQYLRGLVARALAIQNWVDKAEKQALLSETLDLSELFHPDTFLNALRQETARAVGRSVDSLKFVASWKGRLQEAKLQIKISGLLLEGCSFDGNQLSENQLDSPSVSSVLPCFMGWIPQDACGPYSPDECISLPVYTSAERDRVVTNIDVPCGGNQDQWIQCGAALFLKNQ"
    position = 422

    from hommod.tasks import create_model
    path = create_model(sequence, 'HUMAN', position, TemplateID('6F1T', 'e'))

    ok_(path is not None)


@raises(ModelRunError)
@with_setup(setup, end)
def test_create_model_5x0m():
    sequence = "MARFGDEMPARYGGGGSGAAAGVVVGSGGGRGAGGSRQGGQPGAQRMYKQSMAQRARTMALYNPIPVRQNCLTVNRSLFLFSEDNVVRKYAKKITEWPPFEYMILATIIANCIVLALEQHLPDDDKTPMSERLDDTEPYFIGIFCFEAGIKIIALGFAFHKGSYLRNGWNVMDFVVVLTGILATVGTEFDLRTLRAVRVLRPLKLVSGIPSLQVVLKSIMKAMIPLLQIGLLLFFAILIFAIIGLEFYMGKFHTTCFEEGTDDIQGESPAPCGTEEPARTCPNGTKCQPYWEGPNNGITQFDNILFAVLTVFQCITMEGWTDLLYNSNDASGNTWNWLYFIPLIIIGSFFMLNLVLGVLSGEFAKERERVENRRAFLKLRRQQQIERELNGYMEWISKAEEVILAEDETDGEQRHPFDGALRRTTIKKSKTDLLNPEEAEDQLADIASVGSPFARASIKSAKLENSTFFHKKERRMRFYIRRMVKTQAFYWTVLSLVALNTLCVAIVHYNQPEWLSDFLYYAEFIFLGLFMSEMFIKMYGLGTRPYFHSSFNCFDCGVIIGSIFEVIWAVIKPGTSFGISVLRALRLLRIFKVTKYWASLRNLVVSLLNSMKSIISLLFLLFLFIVVFALLGMQLFGGQFNFDEGTPPTNFDTFPAAIMTVFQILTGEDWNEVMYDGIKSQGGVQGGMVFSIYFIVLTLFGNYTLLNVFLAIAVDNLANAQELTKDEQEEEEAANQKLALQKAKEVAEVSPLSAANMSIAVKEQQKNQKPAKSVWEQRTSEMRKQNLLASREALYNEMDPDERWKAAYTRHLRPDMKTHLDRPLVVDPQENRNNNTNKSRAAEPTVDQRLGQQRAEDFLRKQARYHDRARDPSGSAGLDARRPWAGSQEAELSREGPYGRESDHHAREGSLEQPGFWEGEAERGKAGDPHRRHVHRQGGSRESRSGSPRTGADGEHRRHRAHRRPGEEGPEDKAERRARHREGSRPARGGEGEGEGPDGGERRRRHRHGAPATYEGDARREDKERRHRRRKENQGSGVPVSGPNLSTTRPIQQDLGRQDPPLAEDIDNMKNNKLATAESAAPHGSLGHAGLPQSPAKMGNSTDPGPMLAIPAMATNPQNAASRRTPNNPGNPSNPGPPKTPENSLIVTNPSGTQTNSAKTARKPDHTTVDIPPACPPPLNHTVVQVNKNANPDPLPKKEEEKKEEEEDDRGEDGPKPMPPYSSMFILSTTNPLRRLCHYILNLRYFEMCILMVIAMSSIALAAEDPVQPNAPRNNVLRYFDYVFTGVFTFEMVIKMIDLGLVLHQGAYFRDLWNILDFIVVSGALVAFAFTGNSKGKDINTIKSLRVLRVLRPLKTIKRLPKLKAVFDCVVNSLKNVFNILIVYMLFMFIFAVVAVQLFKGKFFHCTDESKEFEKDCRGKYLLYEKNEVKARDREWKKYEFHYDNVLWALLTLFTVSTGEGWPQVLKHSVDATFENQGPSPGYRMEMSIFYVVYFVVFPFFFVNIFVALIIITFQEQGDKMMEEYSLEKNERACIDFAISAKPLTRHMPQNKQSFQYRMWQFVVSPPFEYTIMAMIALNTIVLMMKFYGASVAYENALRVFNIVFTSLFSLECVLKVMAFGILNYFRDAWNIFDFVTVLGSITDILVTEFGNNFINLSFLRLFRAARLIKLLRQGYTIRILLWTFVQSFKALPYVCLLIAMLFFIYAIIGMQVFGNIGIDVEDEDSDEDEFQITEHNNFRTFFQALMLLFRSATGEAWHNIMLSCLSGKPCDKNSGILTRECGNEFAYFYFVSFIFLCSFLMLNLFVAVIMDNFEYLTRDSSILGPHHLDEYVRVWAEYDPAACGRIHYKDMYSLLRVISPPLGLGKKCPHRVACKRLLRMDLPVADDNTVHFNSTLMALIRTALDIKIAKGGADKQQMDAELRKEMMAIWPNLSQKTLDLLVTPHKSTDLTVGKIYAAMMIMEYYRQSKAKKLQAMREEQDRTPLMFQRMEPPSPTQEGGPGQNALPSTQLDPGGALMAHESGLKESPSWVTQRAQEMFQKTGTWSPEQGPPTDMPNSQPNSQSVEMREMGRDGYSDSEHYLPMEGQGRAASMPRLPAENQRRRGRPRGNNLSTISDTSPMKRSASVLGPKARRLDDYSLERVPPEENQRHHQRRRDRSHRASERSLGRYTDVDTGLGTDLSMTTQSGDLPSKERDQERGRPKDRKHRQHHHHHHHHHHPPPPDKDRYAQERPDHGRARARDQRWSRSPSEGREHMAHRQ"

    position = 369

    from hommod.tasks import create_model

    path = create_model(sequence, 'HUMAN', position, TemplateID('5x0m', 'A'))


@with_setup(setup, end)
def test_create_model_1RO6():
    sequence = "MLHVNDLPPPRRHSWICFDVENGPSPGRSPLDPQAGSSSGLVLHAAFPGHSQRRESFLYRSDSDYDLSPKAMSRNSSLPSEQHGDDLIVTPFAQVLASLRSVRNNFTLLTNLHGAPNKRSPAASQAPVSRVSLQEESYQKLAMETLEELDWCLDQLETIQTYRSVSEMASNKFKRMLNRELTHLSEMSRSGNQVSEYISNTFLDKQNDVEIPSPTQKDREKKKKQQLMTQISGVKKLMHSSSLNNTSISRFGVNTENEDHLAKELEDLNKWGLNIFNVAGYSHNRPLTCIMYAIFQERDLLKTFKISSDTFVTYMMTLEDHYHSDVAYHNSLHAADVAQSTHVLLSTPALDAVFTDLEILAAIFAAAIHDVDHPGVSNQFLINTNSELALMYNDESVLENHHLAVGFKLLQEEHCDIFQNLTKKQRQTLRKMVIDMVLATDMSKHMSLLADLKTMVETKKVTSSGVLLLDNYTDRIQVLRNMVHCADLSNPTKSLELYRQWTDRIMEEFFQQGDKERERGMEISPMCDKHTASVEKSQVGFIDYIVHPLWETWADLVQPDAQDILDTLEDNRNWYQSMIPQSPSPPLDERSRDCQGLMEKFQFELTLEEEDSEGPEKEGEGHSYFSSTKTLCVIDPENRDSLEETDIDIATEDKSPIDT"

    position = 254

    from hommod.tasks import create_model

    path = create_model(sequence, 'HUMAN', position, TemplateID('1RO6', 'A'))

    ok_(path is not None)


@with_setup(setup, end)
def test_create_model_6H0A():
    sequence = "MAKLIALTLLGMGLALFRNHQSSYQTRLNALREVQPVELPNCNLVKGIETGSEDLEILPNGLAFISSGLKYPGIKSFNPNSPGKILLMDLNEEDPTVLELGITGSKFDVSSFNPHGISTFTDEDNAMYLLVVNHPDAKSTVELFKFQEEEKSLLHLKTIRHKLLPNLNDIVAVGPEHFYGTNDHYFLDPYLQSWEMYLGLAWSYVVYYSPSEVRVVAEGFDFANGINISPDGKYVYIAELLAHKIHVYEKHANWTLTPLKSLDFNTLVDNISVDPETGDLWVGCHPNGMKIFFYDSENPPASEVLRIQNILTEEPKVTQVYAENGTVLQGSTVASVYKGKLLIGTVFHKALYCEL"

    position = 224

    from hommod.tasks import create_model

    path = create_model(sequence, 'HUMAN', position, TemplateID('6H0A', 'A'))

    ok_(path is not None)

@with_setup(setup, end)
def test_create_model_3srg():
    sequence = "MAKLIALTLLGMGLALFRNHQSSYQTRLNALREVQPVELPNCNLVKGIETGSEDLEILPNGLAFISSGLKYPGIKSFNPNSPGKILLMDLNEEDPTVLELGITGSKFDVSSFNPHGISTFTDEDNAMYLLVVNHPDAKSTVELFKFQEEEKSLLHLKTIRHKLLPNLNDIVAVGPEHFYGTNDHYFLDPYLQSWEMYLGLAWSYVVYYSPSEVRVVAEGFDFANGINISPDGKYVYIAELLAHKIHVYEKHANWTLTPLKSLDFNTLVDNISVDPETGDLWVGCHPNGMKIFFYDSENPPASEVLRIQNILTEEPKVTQVYAENGTVLQGSTVASVYKGKLLIGTVFHKALYCEL"

    position = 224

    from hommod.tasks import create_model

    path = create_model(sequence, 'HUMAN', position, TemplateID('3SRG', 'A'))

    ok_(path is not None)

@with_setup(setup, end)
def test_create_model_lysc():
    sequence = "MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVACAKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV"

    position = 54

    from hommod.tasks import create_model

    path = create_model(sequence, 'HUMAN', position, TemplateID('1C7P', 'A'))

    ok_(path is not None)


@raises(ModelRunError)
@with_setup(setup, end)
def test_create_model_5oyd():
    sequence = "MGPWGWKLRWTVALLLAAAGTAVGDRCERNEFQCQDGKCISYKWVCDGSAECQDGSDESQETCSPKTCSQDEFRCHDGKCISRQFVCDSDRDCLDGSDEASCPVLTCGPASFQCNSSTCIPQLWACDNDPDCEDGSDEWPQRCRGLYVFQGDSSPCSAFEFHCLSGECIHSSWRCDGGPDCKDKSDEENCAVATCRPDEFQCSDGNCIHGSRQCDREYDCKDMSDEVGCVNVTLCEGPNKFKCHSGECITLDKVCNMARDCRDWSDEPIKECGTNECLDNNGGCSHVCNDLKIGYECLCPDGFQLVAQRRCEDIDECQDPDTCSQLCVNLEGGYKCQCEEGFQLDPHTKACKAVGSIAYLFFTNRHEVRKMTLDRSEYTSLIPNLRNVVALDTEVASNRIYWSDLSQRMICSTQLDRAHGVSSYDTVISRDIQAPDGLAVDWIHSNIYWTDSVLGTVSVADTKGVKRKTLFRENGSKPRAIVVDPVHGFMYWTDWGTPAKIKKGGLNGVDIYSLVTENIQWPNGITLDLLSGRLYWVDSKLHSISSIDVNGGNRKTILEDEKRLAHPFSLAVFEDKVFWTDIINEAIFSANRLTGSDVNLLAENLLSPEDMVLFHNLTQPRGVNWCERTTLSNGGCQYLCLPAPQINPHSPKFTCACPDGMLLARDMRSCLTEAEAAVATQETSTVRLKVSSTAVRTQHTTTRPVPDTSRLPGATPGLTTVEIVTMSHQALGDVAGRGNEKKPSSVRALSIVLPIVLLVFLCLGVFLLWKNWRLKNINSINFDNPVYQKTTEDEVHICHNQDGYSYPSRQMVSLEDDVA"

    position = 38

    from hommod.tasks import create_model

    path = create_model(sequence, 'HUMAN', position, TemplateID('5OY9', 'D'))

    ok_(path is None)
