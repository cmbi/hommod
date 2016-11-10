from hommod_rest.factory import create_app, create_celery_app

app = create_app()
celery = create_celery_app(app)
