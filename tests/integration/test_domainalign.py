from hommod_rest.services.domainalign import getAlignments
from hommod_rest.services.modelutils import YasaraChain
from nose.tools import ok_, eq_, with_setup

from hommod_rest.services.interpro import interpro
from hommod_rest.services.blast import blaster
from hommod_rest.services.secstr import secstr
from hommod_rest.services.align import aligner
from hommod_rest.services.model import modeler

import hommod_rest.default_settings as config

import os

mydir = os.path.dirname(__file__)


def setup():
    interpro.interpro_exe = config.INTERPROSCAN
    interpro.storage_dir = config.INTERPRODIR
    blaster.uniprot_db = os.path.join(mydir, "data/mini")  # small version
    blaster.templates_db = config.TEMPLATESDB
    blaster.models_db = config.MODELSDB
    blaster.blastp_exe = config.BLASTP
    secstr.dssp_dir = config.DSSPDIR
    secstr.yasara_dir = config.YASARADIR
    aligner.clustal_exe = config.CLUSTAL
    aligner.kmad_exe = config.KMAD
    modeler.yasara_dir = config.YASARADIR
    modeler.execution_root_dir = config.EXECUTIONDIR
    modeler.model_root_dir = config.MODELDIR
    modeler.template_blacklist = config.TEMPLATE_BLACKLIST


def teardown():
    pass


@with_setup(setup, teardown)
def test_1CBO_A():
    seq = "GYVPAVVIGTGYGAAVSALRLGEAGVQTLMLEMGQLWNQPGPDGNIFCGMLNPDKRSS" \
        + "WFKNRTEAPLGSFLWLDVVNRNIDPYAGVLDRVNYDQMSVYVGRGVGGGSLVNGGMAV" \
        + "EPKRSYFEEILPRVDSSEMYDRYFPRANSMLRVNHIDTKWFEDTEWYKFARVSREQAG" \
        + "KAGLGTVFVPNVYDFGYMQREAAGEVPKSALATEVIYGNNHGKQSLDKTYLAAALGTG" \
        + "KVTIQTLHQVKTIRQTKDGGYALTVEQKDTDGKLLATKEISCRYLFLGAGSLGSTELL" \
        + "VRARDTGTLPNLNSEVGAGWGPNGNIMTARANHMWNPTGAHQSSIPALGIDAWDNSDS" \
        + "SVFAEIAPMPAGLETWVSLYLAITKNPQRGTFVYDAATDRAKLNWTRDQNAPAVNAAK" \
        + "ALFDRINKANGTIYRYDLFGTQLKAFADDFCYNPLGGCVLGKATDDYGRVAGYKNLYV" \
        + "TDGSLIPGSVGVNPFVTITALAERNVERIIKQDV"

    domainranges = interpro.get_domain_locations(seq)
    ok_(len(domainranges) > 0)

    alignments = getAlignments(domainranges, seq)

    eq_(type(alignments), list)
    ok_(len(alignments) > 0)

@with_setup(setup, teardown)
def test_4RH7_A():
    seq = ("MANGTADVRKLFIFTTTQNYFGLMSELWDQPLLCNCLEINNFLDDGNQMLLRVQRSDAGISFSN" +
           "TIEFGDTKDKVLVFFKLRPEVITDENLHDNILVSSMLESPISSLYQAVRQVFAPMLLKDQEWSR" +
           "NFDPKLQNLLSELEAGLGIVLRRSDTNLTKLKFKEDDTRGILTPSDEFQFWIEQAHRGNKQISK" +
           "ERANYFKELFETIAREFYNLDSLSLLEVVDLVETTQDVVDDVWRQTEHDHYPESRMLHLLDIIG" +
           "GSFGRFVQKKLGTLNLWEDPYYLVKESLKAGISICEQWVIVCNHLTGQVWQRYVPHPWKNEKYF" +
           "PETLDKLGKRLEEVLAIRTIHEKFLYFLPASEEKIICLTRVFEPFTGLNPVQYNPYTEPLWKAA" +
           "VSQYEKIIAPAEQKIAGKLKNYISEIQDSPQQLLQAFLKYKELVKRPTISKELMLERETLLARL" +
           "VDSIKDFRLDFENRCRGIPGDASGPLSGKNLSEVVNSIVWVRQLELKVDDTIKIAEALLSDLPG" +
           "FRCFHQSAKDLLDQLKLYEQEQFDDWSRDIQSGLSDSRSGLCIEASSRIMELDSNDGLLKVHYS" +
           "DRLVILLREVRQLSALGFVIPAKIQQVANIAQKFCKQAIILKQVAHFYNSIDQQMIQSQRPMML" +
           "QSALAFEQIIKNSKAGSGGKSQITWDNPKELEGYIQKLQNAAERLATENRKLRKWHTTFCEKVV" +
           "VLMNIDLLRQQQRWKDGLQELRTGLATVEAQGFQASDMHAWKQHWNHQLYKALEHQYQMGLEAL" +
           "NENLPEINIDLTYKQGRLQFRPPFEEIRAKYYREMKRFIGIPNQFKGVGEAGDESIFSIMIDRN" +
           "ASGFLTIFSKAEDLFRRLSAVLHQHKEWIVIGQVDMEALVEKHLFTVHDWEKNFKALKIKGKEV" +
           "ERLPSAVKVDCLNINCNPVKTVIDDLIQKLFDLLVLSLKKSIQAHLHEIDTFVTEAMEVLTIMP" +
           "QSVEEIGDANLQYSKLQERKPEILPLFQEAEDKNRLLRTVAGGGLETISNLKAKWDKFELMMES" +
           "HQLMIKDQIEVMKGNVKSRLQIYYQELEKFKARWDQLKPGDDVIETGQHNTLDKSAKLIKEKKI" +
           "EFDDLEVTRKKLVDDCHHFRLEEPNFSLASSISKDIESCAQIWAFYEEFQQGFQEMANEDWITF" +
           "RTKTYLFEEFLMNWHDRLRKVEEHSVMTVKLQSEVDKYKIVIPILKYVRGEHLSPDHWLDLFRL" +
           "LGLPRGTSLEKLLFGDLLRVADTIVAKAADLKDLNSRAQGEVTIREALRELDLWGVGAVFTLID" +
           "YEDSQSRTMKLIKDWKDIVNQVGDNRCLLQSLKDSPYYKGFEDKVSIWERKLAELDEYLQNLNH" +
           "IQRKWVYLEPIFGRGALPKEQTRFNRVDEDFRSIMTDIKKDNRVTTLTTHAGIRNSLLTILDQL" +
           "QRCQKSLNEFLEEKRSAFPRFYFIGDDDLLEILGQSTNPSVIQSHLKKLFAGINSVCFDEKSKH" +
           "ITAMKSLEGEVVPFKNKVPLSNNVETWLNDLALEMKKTLEQLLKECVTTGRSSQGAVDPSLFPS" +
           "QILCLAEQIKFTEDVENAIKDHSLHQIETQLVNKLEQYTNIDTSSEDPGNTESGILELKLKALI" +
           "LDIIHNIDVVKQLNQIQVHTTEDWAWKKQLRFYMKSDHTCCVQMVDSEFQYTYEYQGNASKLVY" +
           "TPLTDKCYLTLTQAMKMGLGGNPYGPAGTGKTESVKALGGLLGRQVLVFNCDEGIDVKSMGRIF" +
           "VGLVKCGAWGCFDEFNRLEESVLSAVSMQIQTIQDALKNHRTVCELLGKEVEVNSNSGIFITMN" +
           "PAGKGYGGRQKLPDNLKQLFRPVAMSHPDNELIAEVILYSEGFKDAKVLSRKLVAIFNLSRELL" +
           "TPQQHYDWGLRALKTVLRGSGNLLRQLNKSGTTQNANESHIVVQALRLNTMSKFTFTDCTRFDA" +
           "LIKDVFPGIELKEVEYDELSAALKQVFEEANYEIIPNQIKKALELYEQLCQRMGVVIVGPSGAG" +
           "KSTLWRMLRAALCKTGKVVKQYTMNPKAMPRYQLLGHIDMDTREWSDGVLTNSARQVVREPQDV" +
           "SSWIICDGDIDPEWIESLNSVLDDNRLLTMPSGERIQFGPNVNFVFETHDLSCASPATISRMGM" +
           "IFLSDEETDLNSLIKSWLRNQPAEYRNNLENWIGDYFEKALQWVLKQNDYVVETSLVGTVMNGL" +
           "SHLHGCRDHDEFIINLIRGLGGNLNMKSRLEFTKEVFHWARESPPDFHKPMDTYYDSTRGRLAT" +
           "YVLKKPEDLTADDFSNGLTLPVIQTPDMQRGLDYFKPWLSSDTKQPFILVGPEGCGKGMLLRYA" +
           "FSQLRSTQIATVHCSAQTTSRHLLQKLSQTCMVISTNTGRVYRPKDCERLVLYLKDINLPKLDK" +
           "WGTSTLVAFLQQVLTYQGFYDENLEWVGLENIQIVASMSAGGRLGRHKLTTRFTSIVRLCSIDY" +
           "PEREQLQTIYGAYLEPVLHKNLKNHSIWGSSSKIYLLAGSMVQVYEQVRAKFTVDDYSHYFFTP" +
           "CILTQWVLGLFRYDLEGGSSNHPLDYVLEIVAYEARRLFRDKIVGAKELHLFDIILTSVFQGDW" +
           "GSDILDNMSDSFYVTWGARHNSGARAAPGQPLPPHGKPLGKLNSTDLKDVIKKGLIHYGRDNQN" +
           "LDILLFHEVLEYMSRIDRVLSFPGGSLLLAGRSGVGRRTITSLVSHMHGAVLFSPKISRGYELK" +
           "QFKNDLKHVLQLAGIEAQQVVLLLEDYQFVHPTFLEMINSLLSSGEVPGLYTLEELEPLLLPLK" +
           "DQASQDGFFGPVFNYFTYRIQQNLHIVLIMDSANSNFMINCESNPALHKKCQVLWMEGWSNSSM" +
           "KKIPEMLFSETGGGEKYNDKKRKEEKKKNSVDPDFLKSFLLIHESCKAYGATPSRYMTFLHVYS" +
           "AISSSKKKELLKRQSHLQAGVSKLNEAKALVDELNRKAGEQSVLLKTKQDEADAALQMITVSMQ" +
           "DASEQKTELERLKHRIAEEVVKIEERKNKIDDELKEVQPLVNEAKLAVGNIKPESLSEIRSLRM" +
           "PPDVIRDILEGVLRLMGIFDTSWVSMKSFLAKRGVREDIATFDARNISKEIRESVEELLFKNKG" +
           "SFDPKNAKRASTAAAPLAAWVKANIQYSHVLERIHPLETEQAGLESNLKKTEDRKRKLEELLNS" +
           "VGQKVSELKEKFQSRTSEAAKLEAEVSKAQETIKAAEVLINQLDREHKRWNAQVVEITEELATL" +
           "PKRAQLAAAFITYLSAAPESLRKTCLEEWTKSAGLEKFDLRRFLCTESEQLIWKSEGLPSDDLS" +
           "IENALVILQSRVCPFLIDPSSQATEWLKTHLKDSRLEVINQQDSNFITALELAVRFGKTLIIQE" +
           "MDGVEPVLYPLLRRDLVAQGPRYVVQIGDKIIDYNEEFRLFLSTRNPNPFIPPDAASIVTEVNF" +
           "TTTRSGLRGQLLALTIQHEKPDLEEQKTKLLQQEEDKKIQLAKLEESLLETLATSQGNILENKD" +
           "LIESLNQTKASSALIQESLKESYKLQISLDQERDAYLPLAESASKMYFIISDLSKINNMYRFSL" +
           "AAFLRLFQRALQNKQDSENTEQRIQSLISSLQHMVYEYICRCLFKADQLMFALHFVRGMHPELF" +
           "QENEWDTFTGVVVGDMLRKADSQQKIRDQLPSWIDQERSWAVATLKIALPSLYQTLCFEDAALW" +
           "RTYYNNSMCEQEFPSILAKKVSLFQQILVVQALRPDRLQSAMALFACKTLGLKEVSPLPLNLKR" +
           "LYKETLEIEPILIIISPGADPSQELQELANAERSGECYHQVAMGQGQADLAIQMLKECARNGDW" +
           "LCLKNLHLVVSWLPVLEKELNTLQPKDTFRLWLTAEVHPNFTPILLQSSLKITYESPPGLKKNL" +
           "MRTYESWTPEQISKKDNTHRAHALFSLAWFHAACQERRNYIPQGWTKFYEFSLSDLRAGYNIID" +
           "RLFDGAKDVQWEFVHGLLENAIYGGRIDNYFDLRVLQSYLKQFFNSSVIDVFNQRNKKSIFPYS" +
           "VSLPQSCSILDYRAVIEKIPEDDKPSFFGLPANIARSSQRMISSQVISQLRILGRSITAGSKFD" +
           "REIWSNELSPVLNLWKKLNQNSNLIHQKVPPPNDRQGSPILSFIILEQFNAIRLVQSVHQSLAA" +
           "LSKVIRGTTLLSSEVQKLASALLNQKCPLAWQSKWEGPEDPLQYLRGLVARALAIQNWVDKAEK" +
           "QALLSETLDLSELFHPDTFLNALRQETARAVGRSVDSLKFVASWKGRLQEAKLQIKISGLLLEG" +
           "CSFDGNQLSENQLDSPSVSSVLPCFMGWIPQDACGPYSPDECISLPVYTSAERDRVVTNIDVPC" +
           "GGNQDQWIQCGAALFLKNQ")

    tempobj, oligomerisation = modeler._set_template('4RH7')
    yasara_chain = YasaraChain(modeler.yasara, tempobj, 'A')

    domainranges = interpro.get_domain_locations(seq)
    ok_(len(domainranges) > 0)

    alignments = getAlignments(domainranges, seq, yasara_chain)

    covering_alignments = []
    for ran, template_id, alignment in alignments:
        if (ran.start + 1) <= 3373 and (ran.end + 1) >= 3373:
            covering_alignments.append((ran, template_id, alignment))

    ok_(len(covering_alignments) > 0)


@with_setup(setup, teardown)
def test_147b60cfe16160932b7f390a5b952991_3LY6_A():
    seq = ("MGQGEPSQRSTGLAGLYAAPAASPVFIKGSGMDALGIKSCDFQAARNNEEHHTKALSSRRLFVR" +
           "RGQPFTIILYFRAPVRAFLPALKKVALTAQTGEQPSKINRTQATFPISSLGDRKWWSAVVEERD" +
           "AQSWTISVTTPADAVIGHYSLLLQVSGRKQLLLGQFTLLFNPWNREDAVFLKNEAQRMEYLLNQ" +
           "NGLIYLGTADCIQAESWDFGQFEGDVIDLSLRLLSKDKQVEKWSQPVHVARVLGALLHFLKEQR" +
           "VLPTPQTQATQEGALLNKRRGSVPILRQWLTGRGRPVYDGQAWVLAAVACTVLRCLGIPARVVT" +
           "TFASAQGTGGRLLIDEYYNEEGLQNGEGQRGRIWIFQTSTECWMTRPALPQGYDGWQILHPSAP" +
           "NGGGVLGSCDLVVRAVKEGTLGLTPAVSDLFAAINASCVVWKCCEDGTLELTDSNTKYVGNNIS" +
           "TKGVGSDRCEDITQNYKYPEGSLQEKEVLERVEKEKMEREKDNGIRPPSLETASPLYLLLKAPS" +
           "SLPLRGDAQISVTLVNHSEQEKAVQLAIGVQAVHYNGVLAAKLWRKKLHLTLSANLEKIITIGL" +
           "FFSNFERNPPENTFLRLTAMATHSESNLSCFAQEDIAICRPHLAIKMPEKAEQYQPLTASVSLQ" +
           "NSLDAPMEDCVISILGRGLIHRERSYRFRSVWPENTMCAKFQFTPTHVGLQRLTVEVDCNMFQN" +
           "LTNYKSVTVVAPELSA")

    tempobj, oligomerisation = modeler._set_template('3LY6')
    yasara_chain = YasaraChain(modeler.yasara, tempobj, 'A')

    domainranges = interpro.get_domain_locations(seq)
    ok_(len(domainranges) > 0)

    alignments = getAlignments(domainranges, seq, yasara_chain)
    ok_(len(alignments) > 0)

    for ran, template_id, alignment in alignments:
        ok_(ran.end <= len(seq))


@with_setup(setup, teardown)
def test_219431bedd08494fe25aa71ab19af72d_2YPD_A():
    seq = ("MADAAASPVGKRLLLLFADTAASASASAPAAAAASGDPGPALRTRAWRAGTVRAMSGAVPQDLA" +
           "IFVEFDGCNWKQHSWVKVHAEEVIVLLLEGSLVWAPREDPVLLQGIRVSIAQWPALTFTPLVDK" +
           "LGLGSVVPVEYLLDRELRFLSDANGLHLFQMGTDSQNQILLEHAALRETVNALISDQKLQEIFS" +
           "RGPYSVQGHRVKIYQPEGEEGWLYGVVSHQDSITRLMEVSVTESGEIKSVDPRLIHVMLMDNST" +
           "PQSEGGTLKAVKSSKGKKKRESIEGKDGRRRKSASDSGCDPASKKLKGDRGEVDSNGSDGGEAS" +
           "RGPWKGGNASGEPGLDQRAKQPPSTFVPQINRNIRFATYTKENGRTLVVQDEPVGGDTPASFTP" +
           "YSTATGQTPLAPEVGGAENKEAGKTLEQVGQGIVASAAVVTTASSTPNTVRISDTGLAAGTVPE" +
           "KQKGSRSQASGENSRNSILASSGFGAPLPSSSQPLTFGSGRSQSNGVLATENKPLGFSFGCSSA" +
           "QEAQKDTDLSKNLFFQCMSQTLPTSNYFTTVSESLADDSSSRDSFKQSLESLSSGLCKGRSVLG" +
           "TDTKPGSKAGSSVDRKVPAESMPTLTPAFPRSLLNARTPENHENLFLQPPKLSREEPSNPFLAF" +
           "VEKVEHSPFSSFASQASGSSSSATTVTSKVAPSWPESHSSADSASLAKKKPLFITTDSSKLVSG" +
           "VLGSALTSGGPSLSAMGNGRSSSPTSSLTQPIEMPTLSSSPTEERPTVGPGQQDNPLLKTFSNV" +
           "FGRHSGGFLSSPADFSQENKAPFEAVKRFSLDERSLACRQDSDSSTNSDLSDLSDSEEQLQAKT" +
           "GLKGIPEHLMGKLGPNGERSAELLLGKSKGKQAPKGRPRTAPLKVGQSVLKDVSKVKKLKQSGE" +
           "PFLQDGSCINVAPHLHKCRECRLERYRKFKEQEQDDSTVACRFFHFRRLIFTRKGVLRVEGFLS" +
           "PQQSDPDAMNLWIPSSSLAEGIDLETSKYILANVGDQFCQLVMSEKEAMMMVEPHQKVAWKRAV" +
           "RGVREMCDVCETTLFNIHWVCRKCGFGVCLDCYRLRKSRPRSETEEMGDEEVFSWLKCAKGQSH" +
           "EPENLMPTQIIPGTALYNIGDMVHAARGKWGIKANCPCISRQNKSVLRPAVTNGMSQLPSINPS" +
           "ASSGNETTFSGGGGPAPVTTPEPDHVPKADSTDIRSEEPLKTDSSASNSNSELKAIRPPCPDTA" +
           "PPSSALHWLADLATQKAKEETKEAGSLRSVLNKESHSPFGLDSFNSTAKVSPLTPKLFNSLLLG" +
           "PTASNNKTEGSSLRDLLHSGPGKLPQTPLDTGIPFPPVFSTSSAGVKSKASLPNFLDHIIASVV" +
           "ENKKTSDASKRACNLTDTQKEVKEMVMGLNVLDPHTSHSWLCDGRLLCLHDPSNKNNWKIFREC" +
           "WKQGQPVLVSGVHKKLKSELWKPEAFSQEFGDQDVDLVNCRNCAIISDVKVRDFWDGFEIICKR" +
           "LRSEDGQPMVLKLKDWPPGEDFRDMMPTRFEDLMENLPLPEYTKRDGRLNLASRLPSYFVRPDL" +
           "GPKMYNAYGLITAEDRRVGTTNLHLDVSDAVNVMVYVGIPIGEGAHDEEVLKTIDEGDADEVTK" +
           "QRIHDGKEKPGALWHIYAAKDAEKIRELLRKVGEEQGQENPPDHDPIHDQSWYLDQTLRKRLYE" +
           "EYGVQGWAIVQFLGDAVFIPAGAPHQVHNLYSCIKVAEDFVSPEHVKHCFRLTQEFRHLSNTHT" +
           "NHEDKLQVKNIIYHAVKDAVGTLKAHESKLARS")

    tempobj, oligomerisation = modeler._set_template('2YPD')
    yasara_chain = YasaraChain(modeler.yasara, tempobj, 'A')

    domainranges = interpro.get_domain_locations(seq)
    ok_(len(domainranges) > 0)

    alignments = getAlignments(domainranges, seq, yasara_chain)
    ok_(len(alignments) > 0)

    for ran, template_id, alignment in alignments:
        print ran.start, ran.end
        ok_(ran.end <= len(seq))
