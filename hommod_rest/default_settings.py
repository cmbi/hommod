from kombu import Exchange, Queue

DEBUG = False
TESTING = False

# Celery
CELERY_ACCEPT_CONTENT = ['pickle', 'json', 'msgpack', 'yaml']
CELERY_BROKER_URL = 'amqp://guest@localhost'
CELERY_DEFAULT_QUEUE = 'hommod'
CELERYD_CONCURRENCY = 20
CELERY_QUEUES = (
    Queue('hommod', Exchange('hommod'), routing_key='hommod'),
)
CELERY_TRACK_STARTED = True
CELERY_RESULT_BACKEND = 'redis://localhost/1'


# Time it takes for a model to get outdated:
MAX_MODEL_DAYS = 100


# Email logging settings
MAIL_SERVER = "131.174.165.22"
MAIL_SMTP_PORT = 25
MAIL_FROM = "hommod-rest@cmbi.umcn.nl"
MAIL_TO = ["Coos.Baakman@radboudumc.nl", "Jon.Black@radboudumc.nl"]

TEMPLATE_BLACKLIST = "/data/blacklisted_templates"

# Directory where dssp files are stored:
DSSPDIR = '/mnt/cmbi4/dssp/'

# Output directory to store models in:
MODELDIR = '/data/models/'

# Working directory for yasara:
EXECUTIONDIR = '/data/tmp/'

# Directory where interpro files are stored:
INTERPRODIR = '/data/interpro/'

# Executables:
MSA = '/data/prog/pairwiseUBUNTU/MSA' # made by Joanna Lange
BLASTP = '/usr/bin/blastp' # ncbi
CLUSTAL = '/usr/bin/clustalw'
INTERPROSCAN = '/data/prog/interproscan-5.8-49.0/interproscan.sh'

# Yasara installation:
YASARADIR = '/data/prog/Yasara/yasara/'

# Blast databases:
SPECIESDBDIR = '/data/blast/uniprot-species/'
TEMPLATESDB = '/data/blast/templates'


TEMPLATESFASTA = '/data/fasta/templates.fa'
