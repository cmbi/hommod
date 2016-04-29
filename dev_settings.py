import os

DEBUG = True
RETRY_FAILURE = True

CELERY_BROKER_URL = "amqp://guest@amq"
CELERY_RESULT_BACKEND = 'redis://redis/1'

# Paths
MODELDIR = os.path.abspath('./models/')
INTERPRODIR = os.path.abspath('./interpro/')
EXECUTIONDIR = os.path.abspath('./tmp/')

TEMPLATE_BLACKLIST = "/data/blacklisted_templates"

# Directory where dssp files are stored:
DSSPDIR = '/mnt/cmbi4/dssp/'

# Executables:
KMAD = '/deps/hommod-kmad/hommod_kmad' # made by Joanna Lange
BLASTP = '/usr/bin/blastp' # ncbi
CLUSTAL = '/usr/bin/clustalw'
INTERPROSCAN = '/deps/interproscan/interproscan-5.17-56.0/interproscan.sh'

# Yasara installation:
YASARADIR = '/deps/yasara/yasara'

# Blast databases:
UNIPROTDB = '/data/blast/uniprot'
TEMPLATESDB = '/data/blast/templates'


TEMPLATESFASTA = '/data/fasta/templates.fa'
