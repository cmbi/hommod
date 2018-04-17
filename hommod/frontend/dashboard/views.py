import logging
import traceback

from flask import Blueprint, render_template


_log = logging.getLogger(__name__)

bp = Blueprint('dashboard', __name__)


@bp.route('/', methods=['get'])
def index():
    return render_template('index.html')

@bp.route('/model_info/<model_id>/', methods=['get'])
def model_info(model_id):
    return render_template('model.html', model_id=model_id)

@bp.errorhandler(Exception)
def exception_error_handler(error):  # pragma: no cover
    _log.error("Unhandled exception:\n{}".format(traceback.format_exc(error)))
