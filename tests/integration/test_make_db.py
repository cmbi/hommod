from makeXrayPDBFinder2Fasta import gather
from nose.tools import ok_


def test_templates_gather():
    seqs = gather(['1cul'])
    ok_('pdb|1CRN|A' in seqs)
    ok_('pdb|3HQV|A' not in seqs)
    ok_('pdb|1CUL|A' not in seqs)
