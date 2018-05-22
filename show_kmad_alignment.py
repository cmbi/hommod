
from hommod.default_settings import KMAD_EXE, CLUSTALW_EXE
from hommod.controllers.kmad import kmad_aligner
from hommod.controllers.clustal import clustal_aligner


kmad_aligner.kmad_exe = KMAD_EXE
clustal_aligner.clustalw_exe = CLUSTALW_EXE


def gap_equally(ref, ins):
    s = ""
    j = 0
    for i in range(len(ref)):
        if ref[i].isalpha():
            s += ins[j]
            j += 1
        else:
            s += ref[i]
    return s


sequence = "CWAVAVAVGNDGAVAVAVWC"
secstr   = "EEEEEEEE    EEEEEEEE"
target   = "CWAVAVAVAVAVGGGGGGVAVAVAVAVWC"


kmad_alignment = kmad_aligner.align(sequence, secstr, target)
clustal_alignment = clustal_aligner.align({'template': sequence, 'target': target})

print 'kmad'
print kmad_alignment.target_alignment
print gap_equally(kmad_alignment.template_alignment, secstr)
print kmad_alignment.template_alignment

print 'clustal'
print clustal_alignment.aligned_sequences['target']
print gap_equally(clustal_alignment.aligned_sequences['template'], secstr)
print clustal_alignment.aligned_sequences['template']
