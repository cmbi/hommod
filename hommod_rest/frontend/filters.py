import logging
import re

from jinja2 import Markup


_log = logging.getLogger(__name__)

_RE_SUM = re.compile(r"^[\w \.\n']+")
_RE_PAR = re.compile(r"^:param (\w+):([\w ,\.\n']+)", re.M)
_RE_RET = re.compile(r"^:return:([\w ,\.\n']+)", re.M)


def beautify_docstring(docstring):
    result = ""

    # Beautify the summary
    m_sum = _RE_SUM.match(docstring)
    if m_sum:
        result = '{}\n<p>{}</p>'.format(result, m_sum.group(0))

    # Beautify the parameters
    m_par = _RE_PAR.findall(docstring)
    result = '{}\n<p><strong>Parameters:</strong></p>'.format(result)
    result = '{}\n<dl class="dl-horizontal">'.format(result)
    for p in m_par:
        result = '{}\n<dt><code>{}</code></dt><dd>{}</dd>'.format(result,
                                                                  p[0], p[1])
    result = '{}\n</dl>'.format(result)

    # Beautify the return value
    m_ret = _RE_RET.findall(docstring)
    result = '{}\n<p><strong>Returns:</strong></p>'.format(result)
    for rv in m_ret:
        result = '{}\n<p>{}</p>'.format(result, rv)

    _log.debug ("beautify: converted dosctring (length %d) to markup (length %d)"
                % (len (docstring), len (result)))

    return Markup(result)
