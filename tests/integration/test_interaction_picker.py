from hommod_rest.services.interaction import InteractionPicker
from hommod_rest.services.modelutils import YasaraChain
import hommod_rest.default_settings as config

import os
import sys


def test_interaction_pick():

    # The interaction picker is made to accept if one or more of the aligned
    # residues is part of an interaction.

    alignments = {
        'A': {
            'target':
            "MFADRWLFSTNHKDIGTLYLLFGAWAGVLGTALSLLIRAELGQPGNLLGNDHIYNVIVTA" +
            "HAFVMIFFMVMPIMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSLLLLLASAMVEA" +
            "GAGTGWTVYPPLAGNYSHPGASVDLTIFSLHLAGVSSILGAINFITTIINMKPPAMTQYQ" +
            "TPLFVWSVLITAVLLLLSLPVLAAGITMLLTDRNLNTTFFDPAGGGDPILYQHLFWFFGH" +
            "PEVYILILPGFGMISHIVTYYSGKKEPFGYMGMVWAMMSIGFLGFIVWAHHMFTVGMDVD" +
            "TRAYFTSATMIIAIPTGVKVFSWLATLHGSNMKWSAAVLWALGFIFLFTVGGLTGIVLAN" +
            "SSLDIVLHDTYYVVAHFHYVLSMGAVFAIMGGFIHWFPLFSGYTLDQTYAKIHFTIMFIG" +
            "VNLTFFPQHFLGLSGMPRRYSDYPDAYTTWNILSSVGSFISLTAVMLMIFMIWEAFASKR" +
            "KVLMVEEPSMNLEWLYGCPPPYHTFEEPVYM---",
            'template':
            "MFINRWLFSTNHKDIGTLYLLFGAWAGMVGTALSLLIRAELGQPGTLLGDDQIYNVVVTA" +
            "HAFVMIFFMVMPIMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSFLLLLASSMVEA" +
            "GAGTGWTVYPPLAGNLAHAGASVDLTIFSLHLAGVSSILGAINFITTIINMKPPAMSQYQ" +
            "TPLFVWSVMITAVLLLLSLPVLAAGITMLLTDRNLNTTFFDPAGGGDPILYQHLFWFFGH" +
            "PEVYILILPGFGMISHIVTYYSGKKEPFGYMGMVWAMMSIGFLGFIVWAHHMFTVGMDVD" +
            "TRAYFTSATMIIAIPTGVKVFSWLATLHGGNIKWSPAMMWALGFIFLFTVGGLTGIVLAN" +
            "SSLDIVLHDTYYVVAHFHYVLSMGAVFAIMGGFVHWFPLFSGYTLNDTWAKIHFAIMFVG" +
            "VNMTFFPQHFLGLSGMPRRYSDYPDAYTMWNTISSMGSFISLTAVMLMVFIIWEAFASKR" +
            "EVLTVDLTTTNLEWLNGCPPPYHTFEEPTYVNLK"
        },
        # This alignment is OK, the missing helix doesn't interact.
        'B': {
            'target':
            "MAHAAQVGLQDATSPIMEELITFHDHALMIIFLICFLVLYALFLTLTTKLTNTNISDAQE" +
            "METVWTILPAIILVLIALPSLRILY-MTD-E-V--NDPSLTIKSIGHQWYWTYEYTDYGG" +
            "LIFNSYMLPPLFLEPGDLRLLDVDNRVVLPIEAPIRMMITSQDVLHSWAVPTLGLKTDAI" +
            "PGRLNQTTFTATRPGVYYGQCSEICGANHSFMPIVLELI-------------",
            'template':
            "MAYPMQLGFQDATSPIMEELLHFHDHTLMIVFLISSLVLYIISLMLTTKLTHTSTMDAQE" +
            "VETIWTILPAIILILIALPSLRILYMMDDEEIINNNNPSLTVKTMGHQWYWSYEYTDYED" +
            "LSFDSYMIPTSELKPGELRLLEVDNRVVLPMEMTIRMLVSSEDVLHSWAVPSLGLKTDAI" +
            "PGRLNQTTLMSSRPGLYYGQCSEICGSNHSFMPIVLELVPLKYFEKWSASML"
        },
        # This alignment is not OK, the missing helices are the only point of
        # interaction with chain A
        'C': {
            'target':
            "------------------------------------------------------------" +
            "------------------------------------------------------------" +
            "--PLNPLEVPLLNTSVLLASGVSITWAHHSLMENNRNQMIQALLITILLGLYFTLLQASE" +
            "YF---------------------HGLHVIIGSTFLTICFIRQLMFHFTSKHHFGFEAAAW" +
            "YWHFVDVVWLFLYVSIYWWGS",
            'template':
            "MTHQTHAYHMVNPSPWPLTGALSALLMTSGLTMWFHFNSMTLLMIGLTTNMLTMYQWWRD" +
            "VIRESTFQGHHTPAVQKGLRYGMILFIISEVLFFTGFFWAFYHSSLAPTPELGGCWPPTG" +
            "IHPLNPLEVPLLNTSVLLASGVSITWAHHSLMEGDRKHMLQALFITITLGVYFTLLQASE" +
            "YYEAPFTISDGVYGSTFFVATGFHGLHVIIGSTFLIVCFFRQLKFHFTSNHHFGFEAGAW" +
            "YWHFVDVVWLFLYVSIYWWGS"
        },
    }

    sys.path.append(os.path.join(config.YASARADIR, "pym"))
    sys.path.append(os.path.join(config.YASARADIR, "plg"))
    import yasaramodule as yasara
    yasara.info.mode = 'txt'

    obj = yasara.LoadPDB("1OCO", download='yes')[0]

    chains = {}
    for chain_id in ['A', 'B']:
        chains[chain_id] = YasaraChain(yasara, obj, chain_id)

    picker = InteractionPicker('B', chains,
                               {'A': alignments['A']})

    assert picker.accepts('COX2_HUMAN', alignments['B'])

    chains = {}
    for chain_id in ['A', 'C']:
        chains[chain_id] = YasaraChain(yasara, obj, chain_id)

    picker = InteractionPicker('C', chains,
                               {'A': alignments['A']})

    assert (not picker.accepts('COX3_HUMAN', alignments['C']))
