#!/usr/bin/env bash
export HOMMOD_REST_SETTINGS='../dev_settings.py'
gunicorn --log-level debug --log-file "-" -k gevent -b 127.0.0.1:6000 hommod_rest.application:app
