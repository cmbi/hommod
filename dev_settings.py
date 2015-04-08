from kombu import Exchange, Queue
import os

DEBUG = True

# Celery
CELERY_BROKER_URL = 'amqp://guest@localhost'
CELERY_DEFAULT_QUEUE = 'hommod'
CELERY_QUEUES = (
    Queue('hommod', Exchange('hommod'), routing_key='hommod'),
)
CELERY_RESULT_BACKEND = 'amqp'
CELERY_TASK_RESULT_EXPIRES = 18000  # 5 hours.
CELERY_TRACK_STARTED = True

MAX_MODEL_DAYS = 100
RETRY_FAILURE = True

# Paths
MODELDIR = os.path.abspath('./models/')
INTERPRODIR = os.path.abspath('./interpro/')
EXECUTIONDIR = os.path.abspath('./tmp/')

BLASTP = '/usr/local/blastp'
CLUSTAL = '/usr/local/clustalw'
MSA = '/data/prog/pairwiseUBUNTU/MSA'
INTERPROSCAN = '/data/prog/interproscan-5.10-50.0/interproscan.sh'
YASARADIR = '/data/prog/Yasara/yasara/'
SPECIESDBDIR = '/data/blast/uniprot-species/'
TEMPLATESDB = '/data/blast/templates/templates'
DSSPDIR = '/mnt/cmbi4/dssp/'
