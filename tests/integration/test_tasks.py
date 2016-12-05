from hommod_rest.factory import create_app, create_celery_app
from hommod_rest.services.modelutils import TemplateID

import logging
_log = logging.getLogger(__name__)


class TestTasks:

    @classmethod
    def setup_class(self):
        self.app = create_app({'TESTING': True, 'CELERY_ALWAYS_EAGER': True})
        self.celery = create_celery_app(self.app)

    def test_create_model(self):
        # Best pick a simple model which doesn't take much time ..
        from hommod_rest.tasks import create_model
        path = create_model("IICCASITARSDFDVCRLAGTAQAVCAVFVGCVIVPAATCPGDFGD",
                            "CRAAB", 25, None)

        _log.debug("Recieved model path: \'%s\'" % path)

        assert (len(path) > 0)

    def test_create_model_4RH7_A(self):
        # This is a slow model!
        from hommod_rest.tasks import create_model
        path = create_model(
            "MANGTADVRKLFIFTTTQNYFGLMSELWDQPLLCNCLEINNFLDDGNQMLLRVQRSDAGISFSNTIEFGDTK" +
            "DKVLVFFKLRPEVITDENLHDNILVSSMLESPISSLYQAVRQVFAPMLLKDQEWSRNFDPKLQNLLSELEAG" +
            "LGIVLRRSDTNLTKLKFKEDDTRGILTPSDEFQFWIEQAHRGNKQISKERANYFKELFETIAREFYNLDSLS" +
            "LLEVVDLVETTQDVVDDVWRQTEHDHYPESRMLHLLDIIGGSFGRFVQKKLGTLNLWEDPYYLVKESLKAGI" +
            "SICEQWVIVCNHLTGQVWQRYVPHPWKNEKYFPETLDKLGKRLEEVLAIRTIHEKFLYFLPASEEKIICLTR" +
            "VFEPFTGLNPVQYNPYTEPLWKAAVSQYEKIIAPAEQKIAGKLKNYISEIQDSPQQLLQAFLKYKELVKRPT" +
            "ISKELMLERETLLARLVDSIKDFRLDFENRCRGIPGDASGPLSGKNLSEVVNSIVWVRQLELKVDDTIKIAE" +
            "ALLSDLPGFRCFHQSAKDLLDQLKLYEQEQFDDWSRDIQSGLSDSRSGLCIEASSRIMELDSNDGLLKVHYS" +
            "DRLVILLREVRQLSALGFVIPAKIQQVANIAQKFCKQAIILKQVAHFYNSIDQQMIQSQRPMMLQSALAFEQ" +
            "IIKNSKAGSGGKSQITWDNPKELEGYIQKLQNAAERLATENRKLRKWHTTFCEKVVVLMNIDLLRQQQRWKD" +
            "GLQELRTGLATVEAQGFQASDMHAWKQHWNHQLYKALEHQYQMGLEALNENLPEINIDLTYKQGRLQFRPPF" +
            "EEIRAKYYREMKRFIGIPNQFKGVGEAGDESIFSIMIDRNASGFLTIFSKAEDLFRRLSAVLHQHKEWIVIG" +
            "QVDMEALVEKHLFTVHDWEKNFKALKIKGKEVERLPSAVKVDCLNINCNPVKTVIDDLIQKLFDLLVLSLKK" +
            "SIQAHLHEIDTFVTEAMEVLTIMPQSVEEIGDANLQYSKLQERKPEILPLFQEAEDKNRLLRTVAGGGLETI" +
            "SNLKAKWDKFELMMESHQLMIKDQIEVMKGNVKSRLQIYYQELEKFKARWDQLKPGDDVIETGQHNTLDKSA" +
            "KLIKEKKIEFDDLEVTRKKLVDDCHHFRLEEPNFSLASSISKDIESCAQIWAFYEEFQQGFQEMANEDWITF" +
            "RTKTYLFEEFLMNWHDRLRKVEEHSVMTVKLQSEVDKYKIVIPILKYVRGEHLSPDHWLDLFRLLGLPRGTS" +
            "LEKLLFGDLLRVADTIVAKAADLKDLNSRAQGEVTIREALRELDLWGVGAVFTLIDYEDSQSRTMKLIKDWK" +
            "DIVNQVGDNRCLLQSLKDSPYYKGFEDKVSIWERKLAELDEYLQNLNHIQRKWVYLEPIFGRGALPKEQTRF" +
            "NRVDEDFRSIMTDIKKDNRVTTLTTHAGIRNSLLTILDQLQRCQKSLNEFLEEKRSAFPRFYFIGDDDLLEI" +
            "LGQSTNPSVIQSHLKKLFAGINSVCFDEKSKHITAMKSLEGEVVPFKNKVPLSNNVETWLNDLALEMKKTLE" +
            "QLLKECVTTGRSSQGAVDPSLFPSQILCLAEQIKFTEDVENAIKDHSLHQIETQLVNKLEQYTNIDTSSEDP" +
            "GNTESGILELKLKALILDIIHNIDVVKQLNQIQVHTTEDWAWKKQLRFYMKSDHTCCVQMVDSEFQYTYEYQ" +
            "GNASKLVYTPLTDKCYLTLTQAMKMGLGGNPYGPAGTGKTESVKALGGLLGRQVLVFNCDEGIDVKSMGRIF" +
            "VGLVKCGAWGCFDEFNRLEESVLSAVSMQIQTIQDALKNHRTVCELLGKEVEVNSNSGIFITMNPAGKGYGG" +
            "RQKLPDNLKQLFRPVAMSHPDNELIAEVILYSEGFKDAKVLSRKLVAIFNLSRELLTPQQHYDWGLRALKTV" +
            "LRGSGNLLRQLNKSGTTQNANESHIVVQALRLNTMSKFTFTDCTRFDALIKDVFPGIELKEVEYDELSAALK" +
            "QVFEEANYEIIPNQIKKALELYEQLCQRMGVVIVGPSGAGKSTLWRMLRAALCKTGKVVKQYTMNPKAMPRY" +
            "QLLGHIDMDTREWSDGVLTNSARQVVREPQDVSSWIICDGDIDPEWIESLNSVLDDNRLLTMPSGERIQFGP" +
            "NVNFVFETHDLSCASPATISRMGMIFLSDEETDLNSLIKSWLRNQPAEYRNNLENWIGDYFEKALQWVLKQN" +
            "DYVVETSLVGTVMNGLSHLHGCRDHDEFIINLIRGLGGNLNMKSRLEFTKEVFHWARESPPDFHKPMDTYYD" +
            "STRGRLATYVLKKPEDLTADDFSNGLTLPVIQTPDMQRGLDYFKPWLSSDTKQPFILVGPEGCGKGMLLRYA" +
            "FSQLRSTQIATVHCSAQTTSRHLLQKLSQTCMVISTNTGRVYRPKDCERLVLYLKDINLPKLDKWGTSTLVA" +
            "FLQQVLTYQGFYDENLEWVGLENIQIVASMSAGGRLGRHKLTTRFTSIVRLCSIDYPEREQLQTIYGAYLEP" +
            "VLHKNLKNHSIWGSSSKIYLLAGSMVQVYEQVRAKFTVDDYSHYFFTPCILTQWVLGLFRYDLEGGSSNHPL" +
            "DYVLEIVAYEARRLFRDKIVGAKELHLFDIILTSVFQGDWGSDILDNMSDSFYVTWGARHNSGARAAPGQPL" +
            "PPHGKPLGKLNSTDLKDVIKKGLIHYGRDNQNLDILLFHEVLEYMSRIDRVLSFPGGSLLLAGRSGVGRRTI" +
            "TSLVSHMHGAVLFSPKISRGYELKQFKNDLKHVLQLAGIEAQQVVLLLEDYQFVHPTFLEMINSLLSSGEVP" +
            "GLYTLEELEPLLLPLKDQASQDGFFGPVFNYFTYRIQQNLHIVLIMDSANSNFMINCESNPALHKKCQVLWM" +
            "EGWSNSSMKKIPEMLFSETGGGEKYNDKKRKEEKKKNSVDPDFLKSFLLIHESCKAYGATPSRYMTFLHVYS" +
            "AISSSKKKELLKRQSHLQAGVSKLNEAKALVDELNRKAGEQSVLLKTKQDEADAALQMITVSMQDASEQKTE" +
            "LERLKHRIAEEVVKIEERKNKIDDELKEVQPLVNEAKLAVGNIKPESLSEIRSLRMPPDVIRDILEGVLRLM" +
            "GIFDTSWVSMKSFLAKRGVREDIATFDARNISKEIRESVEELLFKNKGSFDPKNAKRASTAAAPLAAWVKAN" +
            "IQYSHVLERIHPLETEQAGLESNLKKTEDRKRKLEELLNSVGQKVSELKEKFQSRTSEAAKLEAEVSKAQET" +
            "IKAAEVLINQLDREHKRWNAQVVEITEELATLPKRAQLAAAFITYLSAAPESLRKTCLEEWTKSAGLEKFDL" +
            "RRFLCTESEQLIWKSEGLPSDDLSIENALVILQSRVCPFLIDPSSQATEWLKTHLKDSRLEVINQQDSNFIT" +
            "ALELAVRFGKTLIIQEMDGVEPVLYPLLRRDLVAQGPRYVVQIGDKIIDYNEEFRLFLSTRNPNPFIPPDAA" +
            "SIVTEVNFTTTRSGLRGQLLALTIQHEKPDLEEQKTKLLQQEEDKKIQLAKLEESLLETLATSQGNILENKD" +
            "LIESLNQTKASSALIQESLKESYKLQISLDQERDAYLPLAESASKMYFIISDLSKINNMYRFSLAAFLRLFQ" +
            "RALQNKQDSENTEQRIQSLISSLQHMVYEYICRCLFKADQLMFALHFVRGMHPELFQENEWDTFTGVVVGDM" +
            "LRKADSQQKIRDQLPSWIDQERSWAVATLKIALPSLYQTLCFEDAALWRTYYNNSMCEQEFPSILAKKVSLF" +
            "QQILVVQALRPDRLQSAMALFACKTLGLKEVSPLPLNLKRLYKETLEIEPILIIISPGADPSQELQELANAE" +
            "RSGECYHQVAMGQGQADLAIQMLKECARNGDWLCLKNLHLVVSWLPVLEKELNTLQPKDTFRLWLTAEVHPN" +
            "FTPILLQSSLKITYESPPGLKKNLMRTYESWTPEQISKKDNTHRAHALFSLAWFHAACQERRNYIPQGWTKF" +
            "YEFSLSDLRAGYNIIDRLFDGAKDVQWEFVHGLLENAIYGGRIDNYFDLRVLQSYLKQFFNSSVIDVFNQRN" +
            "KKSIFPYSVSLPQSCSILDYRAVIEKIPEDDKPSFFGLPANIARSSQRMISSQVISQLRILGRSITAGSKFD" +
            "REIWSNELSPVLNLWKKLNQNSNLIHQKVPPPNDRQGSPILSFIILEQFNAIRLVQSVHQSLAALSKVIRGT" +
            "TLLSSEVQKLASALLNQKCPLAWQSKWEGPEDPLQYLRGLVARALAIQNWVDKAEKQALLSETLDLSELFHP" +
            "DTFLNALRQETARAVGRSVDSLKFVASWKGRLQEAKLQIKISGLLLEGCSFDGNQLSENQLDSPSVSSVLPC" +
            "FMGWIPQDACGPYSPDECISLPVYTSAERDRVVTNIDVPCGGNQDQWIQCGAALFLKNQ",
            "HUMAN", 3373, TemplateID('4RH7', 'A'))

        assert (len(path) > 0)

    def test_take_template(self):
        # Best pick a simple model which doesn't take much time ..
        from hommod_rest.tasks import create_model
        path = create_model("TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN",
                            "CRAAB", 25, None)

        _log.debug("Recieved model path: \'%s\'" % path)

        assert (len(path) > 0)

    def test_model_high_memory(self):
        from hommod_rest.tasks import create_model
        path = create_model(
                "QRNLRKYLQLRTWPWYKLWQKVKPLLNVSRIEDEIARLEEKAKKAEELHAAEVKVRKELE",
                "DROME", 820, TemplateID("3JBH", 'A'))
        assert (len(path) > 0)

    def test_model_1R6F(self):
        from hommod_rest.tasks import create_model
        path = create_model(
            "MERAESSSTEPAKAIKPIDRKSVHQICSGQVVLSLSTAVKELVENSLDAGATNIDLKLKDYG" +
            "VDLIEVSDNGCGVEEENFEGLTLKHHTSKIQEFADLTQVETFGFRGEALSSLCALSDVTIST" +
            "CHASAKVGTRLMFDHNGKIIQKTPYPRPRGTTVSVQQLFSTLPVRHKEFQRNIKKEYAKMVQ" +
            "VLHAYCIISAGIRVSCTNQLGQGKRQPVVCTGGSPSIKENIGSVFGQKQLQSLIPFVQLPPS" +
            "DSVCEEYGLSCSDALHNLFYISGFISQCTHGVGRSSTDRQFFFINRRPCDPAKVCRLVNEVY" +
            "HMYNRHQYPFVVLNISVDSECVDINVTPDKRQILLQEEKLLLAVLKTSLIGMFDSDVNKLNV" +
            "SQQPLLDVEGNLIKMHAADLEKPMVEKQDQSPSLRTGEEKKDVSISRLREAFSLRHTTENKP" +
            "HSPKTPEPRRSPLGQKRGMLSSSTSGAISDKGVLRPQKEAVSSSHGPSDPTDRAEVEKDSGH" +
            "GSTSVDSEGFSIPDTGSHCSSEYAASSPGDRGSQEHVDSQEKAPETDDSFSDVDCHSNQEDT" +
            "GCKFRVLPQPTNLATPNTKRFKKEEILSSSDICQKLVNTQDMSASQVDVAVKINKKVVPLDF" +
            "SMSSLAKRIKQLHHEAQQSEGEQNYRKFRAKICPGENQAAEDELRKEISKTMFAEMEIIGQF" +
            "NLGFIITKLNEDIFIVDQHATDEKYNFEMLQQHTVLQGQRLIAPQTLNLTAVNEAVLIENLE" +
            "IFRKNGFDFVIDENAPVTERAKLISLPTSKNWTFGPQDVDELIFMLSDSPGVMCRPSRVKQM" +
            "FASRACRKSVMIGTALNTSEMKKLITHMGEMDHPWNCPHGRPTMRHIANLGVISQN",
            "HUMAN", 470, TemplateID("1R6F", 'A'))
        # Alignment is poor, shouldn't produce any model!
        assert (len(path) <= 0)

