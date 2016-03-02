import logging
from logging.handlers import SMTPHandler

from celery import Celery
from flask import Flask

_log = logging.getLogger(__name__)


def create_app(settings=None):
    _log.info("Creating app")

    app = Flask(__name__)
    app.config.from_object('hommod_rest.default_settings')
    if settings:
        app.config.update(settings)
    else:  # pragma: no cover
        app.config.from_envvar('HOMMOD_REST_SETTINGS')  # pragma: no cover

    # Ignore Flask's built-in logging
    # app.logger is accessed here so Flask tries to create it
    app.logger_name = "nowhere"
    app.logger

    # Configure logging.
    #
    # It is somewhat dubious to get _log from the root package, but I can't see
    # a better way. Having the email handler configured at the root means all
    # child loggers inherit it.
    from hommod_rest import _log as hommod_logger

    # Only log to email during production.
    if not app.debug and not app.testing:  # pragma: no cover
        mail_handler = SMTPHandler((app.config["MAIL_SERVER"],
                                   app.config["MAIL_SMTP_PORT"]),
                                   app.config["MAIL_FROM"],
                                   app.config["MAIL_TO"],
                                   "hommod-rest failed")
        mail_handler.setLevel(logging.ERROR)
        hommod_logger.addHandler(mail_handler)
        mail_handler.setFormatter(
            logging.Formatter("Message type: %(levelname)s\n" +
                              "Location: %(pathname)s:%(lineno)d\n" +
                              "Module: %(module)s\n" +
                              "Function: %(funcName)s\n" +
                              "Time: %(asctime)s\n" +
                              "Message:\n" +
                              "%(message)s"))

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

    # Use ProxyFix to correct URL's when redirecting.
    from hommod_rest.middleware import ReverseProxied
    app.wsgi_app = ReverseProxied(app.wsgi_app)

    # Register blueprints
    from hommod_rest.frontend.api.endpoints import bp as api_bp
    app.register_blueprint(api_bp)

    # Initialise services
    from hommod_rest.services.model import modeler
    modeler.yasara_dir = app.config ['YASARADIR']
    modeler.execution_root_dir = app.config ['EXECUTIONDIR']
    modeler.model_root_dir = app.config ['MODELDIR']
    modeler.template_blacklist = app.config ['TEMPLATE_BLACKLIST']

    from hommod_rest.services.secstr import secstr
    secstr.dssp_dir = app.config ['DSSPDIR']
    secstr.yasara_dir = app.config ['YASARADIR']

    from hommod_rest.services.interpro import interpro
    interpro.interproExe = app.config ['INTERPROSCAN']
    interpro.storageDir = app.config ['INTERPRODIR']

    from hommod_rest.services.blast import blaster
    blaster.templatesDB = app.config ['TEMPLATESDB']
    blaster.uniprotspeciesDir = app.config ['SPECIESDBDIR']
    blaster.blastpExe = app.config ['BLASTP']

    from hommod_rest.services.align import aligner
    aligner.clustalExe = app.config ['CLUSTAL']
    aligner.msaExe = app.config ['MSA']

    import hommod_rest.services.modelutils
    hommod_rest.services.modelutils.TEMPLATESFASTA = app.config ['TEMPLATESFASTA']

    return app


def create_celery_app(app):  # pragma: no cover
    _log.info("Creating celery app")

    app = app or create_app()

    celery = Celery(__name__, backend=app.config['CELERY_RESULT_BACKEND'],
                    broker=app.config['CELERY_BROKER_URL'])
    celery.conf.update(app.config)
    TaskBase = celery.Task

    class ContextTask(TaskBase):
        abstract = True

        def __call__(self, *args, **kwargs):
            with app.app_context():
                return TaskBase.__call__(self, *args, **kwargs)

    celery.Task = ContextTask

    import hommod_rest.tasks

    return celery
