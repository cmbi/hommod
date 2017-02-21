from kombu import Exchange, Queue
from celery.schedules import crontab

DEBUG = False
TESTING = False

# Celery
CELERY_ACCEPT_CONTENT = ['pickle', 'json', 'msgpack', 'yaml']
CELERY_BROKER_URL = 'amqp://guest@hommodrest_rabbitmq_1'
CELERY_DEFAULT_QUEUE = 'hommod'
CELERYD_CONCURRENCY = 20
CELERY_QUEUES = (
    Queue('hommod', Exchange('hommod'), routing_key='hommod'),
)
CELERY_TRACK_STARTED = True
CELERY_RESULT_BACKEND = 'redis://hommodrest_redis_1/1'
CELERYBEAT_SCHEDULE = {
    # 1 Every Hour
    'remodel_oldest_hg': {
        'task': 'hommod_rest.tasks.remodel_oldest_hg',
        'schedule': crontab(minute=0, hour='*'),
    },
}

# Time it takes for a model to get outdated:
MAX_MODEL_DAYS = 100


# Email logging settings
MAIL_SERVER = "131.174.165.22"
MAIL_SMTP_PORT = 25
MAIL_FROM = "hommod-rest@cmbi.umcn.nl"
MAIL_TO = ["Coos.Baakman@radboudumc.nl", "Jon.Black@radboudumc.nl"]

# Folders
TEMPLATE_BLACKLIST = "/data/blacklisted_templates"
PDBFINDER2 = '/mnt/cmbi4/pdbfinder2/PDBFIND2.TXT'
DSSPDIR = '/mnt/cmbi4/dssp/'
HGFASTADIR = '/mnt/cmbi4/hg-fasta/'
MODELDIR = '/data/models/'
EXECUTIONDIR = '/data/tmp/'  # yasara execution dir
INTERPRODIR = '/data/interpro/'

# Executables
KMAD = '/deps/hommod-kmad/hommod_kmad'  # made by Joanna Lange
BLASTP = '/usr/bin/blastp'  # ncbi
CLUSTAL = '/usr/bin/clustalw'
INTERPROSCAN = '/deps/interproscan/interproscan-5.17-56.0/interproscan.sh'
YASARADIR = '/deps/yasara/yasara/'

# Blast databases
TEMPLATESDB = '/data/blast/templates'
UNIPROTDB = '/data/blast/uniprot'
TEMPLATESFASTA = '/data/fasta/templates.fa'
