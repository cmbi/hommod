import os

DEBUG = True
RETRY_FAILURE = True

CELERY_BROKER_URL = "amqp://guest@rabbitmq_1"

# Paths
MODELDIR = os.path.abspath('./models/')
INTERPRODIR = os.path.abspath('./interpro/')
EXECUTIONDIR = os.path.abspath('./tmp/')

TEMPLATE_BLACKLIST = "/data/blacklisted_templates"

# Directory where dssp files are stored:
DSSPDIR = '/mnt/cmbi4/dssp/'

# Executables:
MSA = 'hommod-kmad/hommod_kmad' # made by Joanna Lange
BLASTP = '/usr/bin/blastp' # ncbi
CLUSTAL = '/usr/bin/clustalw'
INTERPROSCAN = '/data/prog/interproscan-5.8-49.0/interproscan.sh'

# Yasara installation:
YASARADIR = '/deps/yasara/yasara'

# Blast databases:
SPECIESDBDIR = '/data/blast/uniprot-species/'
TEMPLATESDB = '/data/blast/templates'


TEMPLATESFASTA = '/data/fasta/templates.fa'
