from kombu import Exchange, Queue

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

MAX_MODEL_DAYS = 100

# Email logging settings
MAIL_SERVER = "131.174.165.22"
MAIL_SMTP_PORT = 25
MAIL_FROM = "hommod-rest@cmbi.umcn.nl"
MAIL_TO = ["Coos.Baakman@radboudumc.nl"]

# Paths

DSSPDIR = '/mnt/cmbi4/dssp/'
MODELDIR = '/data/models/'
EXECUTIONDIR = '/data/tmp/'
INTERPRODIR = '/data/interpro/'

MSA = '/data/prog/pairwiseUBUNTU/MSA'
BLASTP = '/usr/bin/blastp'
CLUSTAL = '/usr/bin/clustalw'
INTERPROSCAN = '/data/prog/interproscan-5.8-49.0/interproscan.sh'
YASARADIR = '/data/prog/Yasara/yasara/'
SPECIESDBDIR = '/data/blast/uniprot-species/'
TEMPLATESDB = '/data/blast/templates/templates'
