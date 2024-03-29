import os
import logging
from logging.handlers import RotatingFileHandler

from celery import Celery
from flask import Flask


_log = logging.getLogger(__name__)
sh = logging.StreamHandler()
file_handler = RotatingFileHandler(os.environ["LOG_FILENAME"], maxBytes=10485760)


def create_app(settings=None):
    app = Flask(__name__, static_folder='frontend/static',
                template_folder='frontend/templates')

    app.config.from_object('hommod.default_settings')
    if settings:
        app.config.update(settings)
    else:  # pragma: no cover
        app.config.from_envvar('HOMMOD_SETTINGS')  # pragma: no cover

    # Ignore Flask's built-in logging
    # app.logger is accessed here so Flask tries to create it
    app.logger_name = "nowhere"
    app.logger

    # Configure logging.
    #
    # It is somewhat dubious to get _log from the root package, but I can't see
    # a better way. Having the email handler configured at the root means all
    # child loggers inherit it.
    from hommod import _log as hommod_logger

    # Only log to the console during development and production, but not during
    # testing.
    if app.testing:
        hommod_logger.setLevel(logging.DEBUG)
    else:
        ch = logging.StreamHandler()
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        hommod_logger.addHandler(ch)

        if app.debug:
            hommod_logger.setLevel(logging.DEBUG)
        else:
            hommod_logger.setLevel(logging.INFO)

            # Log to file in production as well.  Log filename is loaded from
            # the environment, not from the settings file. maxBytes = 10MB.
            file_handler.setLevel(logging.INFO)
            file_handler.setFormatter(formatter)
            hommod_logger.addHandler(file_handler)

    # Use ProxyFix to correct URL's when redirecting.
    from hommod.middleware import ReverseProxied
    app.wsgi_app = ReverseProxied(app.wsgi_app)

    # Register blueprints
    from hommod.frontend.api.endpoints import bp as api_bp
    app.register_blueprint(api_bp)
    from hommod.frontend.dashboard.views import bp as dashboard_bp
    app.register_blueprint(dashboard_bp)

    # Register jinja2 filters
    from hommod.frontend.filters import beautify_docstring
    app.jinja_env.filters['beautify_docstring'] = beautify_docstring

    return app


def create_celery_app(flask_app=None):  # pragma: no cover

    if flask_app:
        _log.info("Creating celery app with given flask app")
    else:
        _log.info("Creating celery app with new flask app")
        flask_app = create_app()

    from hommod import default_settings

    celery = Celery(__name__)
    celery.config_from_object(default_settings)
    celery.conf.update(flask_app.config)
    TaskBase = celery.Task

    class ContextTask(TaskBase):
        abstract = True

        def __call__(self, *args, **kwargs):
            with flask_app.app_context():
                return TaskBase.__call__(self, *args, **kwargs)

    celery.Task = ContextTask

    import hommod.tasks

    # Apply config settings
    from hommod.controllers.uniprot import uniprot
    uniprot.fasta_paths = [flask_app.config['SPROT_FASTA'], flask_app.config['TREMBL_FASTA']]

    from hommod.services.interpro import interpro
    interpro.url = flask_app.config['INTERPRO_URL']
    interpro.email = flask_app.config['ADMIN_EMAIL']

    from hommod.controllers.kmad import kmad_aligner
    kmad_aligner.kmad_exe = flask_app.config['KMAD_EXE']

    from hommod.controllers.clustal import clustal_aligner
    clustal_aligner.clustalw_exe = flask_app.config['CLUSTALW_EXE']

    from hommod.controllers.storage import model_storage
    model_storage.model_dir = flask_app.config['MODEL_DIR']

    from hommod.controllers.model import modeler
    modeler.yasara_dir = flask_app.config['YASARA_DIR']
    modeler.uniprot_databank = flask_app.config['UNIPROT_BLAST_DATABANK']

    from hommod.controllers.domain import domain_aligner
    domain_aligner.forbidden_interpro_domains = flask_app.config['FORBIDDEN_INTERPRO_DOMAINS']
    domain_aligner.similar_ranges_min_overlap_percentage = flask_app.config['SIMILAR_RANGES_MIN_OVERLAP_PERCENTAGE']
    domain_aligner.similar_ranges_max_length_difference_percentage = flask_app.config['SIMILAR_RANGES_MAX_LENGTH_DIFFERENCE_PERCENTAGE']
    domain_aligner.min_percentage_coverage = flask_app.config['DOMAIN_MIN_PERCENTAGE_COVERAGE']
    domain_aligner.template_blast_databank = flask_app.config['TEMPLATE_BLAST_DATABANK']
    domain_aligner.highly_homologous_percentage_identity = flask_app.config['HIGHLY_HOMOLOGOUS_PERCENTAGE_IDENTITY']

    from hommod.controllers.blast import blaster
    blaster.blastp_exe = flask_app.config['BLASTP_EXE']

    from hommod.controllers.blacklist import blacklister
    blacklister.file_path = flask_app.config['BLACKLIST_FILE_PATH']

    from hommod.services.dssp import dssp
    dssp.dssp_dir = flask_app.config['DSSP_DIR']

    from hommod.services.helpers.cache import cache_manager as cm
    cm.redis_hostname = flask_app.config['CACHE_REDIS_HOST']
    cm.redis_port = flask_app.config['CACHE_REDIS_PORT']
    cm.redis_db = flask_app.config['CACHE_REDIS_DB']
    cm.expiration_time = flask_app.config['CACHE_EXPIRATION_TIME']
    cm.lock_timeout = flask_app.config['CACHE_LOCK_TIMEOUT']

    return celery
