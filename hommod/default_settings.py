from kombu import Exchange, Queue
from celery.schedules import crontab
import datetime

DEBUG = False
TESTING = False

# Cache
CACHE_REDIS_HOST = 'hommod_redis_1'
CACHE_REDIS_PORT = 6379
CACHE_REDIS_DB = 1
CACHE_EXPIRATION_TIME = 60*60*24*30  # 30 days
CACHE_LOCK_TIMEOUT = 60*60  # 1 hour

# Celery
CELERY_ACCEPT_CONTENT = ['pickle', 'json', 'msgpack', 'yaml']
CELERY_BROKER_URL = 'amqp://guest@hommod_rabbitmq_1'
CELERY_DEFAULT_QUEUE = 'hommod'
CELERYD_CONCURRENCY = 20
CELERY_QUEUES = ( 
    Queue('hommod', Exchange('hommod'), routing_key='hommod'),
)
CELERY_TRACK_STARTED = True
CELERY_RESULT_BACKEND = 'redis://hommod_redis_1/1'

# Time it takes for a model to get outdated:
MAX_MODEL_DAYS = 100 

# Services
INTERPRO_URL = 'http://www.ebi.ac.uk/Tools/services/rest/iprscan5'
UNIPROT_URL = 'http://www.uniprot.org/uniprot'

# Directories and File Paths
YASARA_DIR = '/deps/yasara/yasara'
MODEL_DIR = '/data/models/'
BLACKLIST_FILE_PATH = '/data/blacklisted_templates'
DSSP_DIR = '/mnt/chelonium/dssp/'
PDBFINDER2_FILE_PATH = '/mnt/chelonium/pdbfinder2/PDBFIND2.TXT'

# Executables
KMAD_EXE = '/deps/hommod-kmad/hommod_kmad'  # made by Joanna Lange
BLASTP_EXE = '/usr/bin/blastp'  # ncbi
CLUSTALW_EXE = '/usr/bin/clustalw'

# Databanks
TEMPLATE_BLAST_DATABANK = '/data/blast/templates'
UNIPROT_BLAST_DATABANK = '/data/blast/uniprot'

# Fastas
SPROT_FASTA = '/data/fasta/uniprot_sprot.fasta'

# Domain alignment settings
HIGHLY_HOMOLOGOUS_PERCENTAGE_IDENTITY = 80.0
DOMAIN_MAX_MERGE_DISTANCE = 10  # amino acids
DOMAIN_MIN_PERCENTAGE_COVERAGE = 80.0
SIMILAR_RANGES_MIN_OVERLAP_PERCENTAGE = 80.0
SIMILAR_RANGES_MAX_LENGTH_DIFFERENCE_PERCENTAGE = 10.0
FORBIDDEN_INTERPRO_DOMAINS = ['IPR003596']  # Ig variable domain like
