#!/usr/bin/env bash
export HOMMOD_REST_SETTINGS='../dev_settings.py'
celery -A hommod_rest.celery_application:celery worker -B
