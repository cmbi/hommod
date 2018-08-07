from threading import local


class ModelLogger:
    """Class for storing information that needs to be writter to the job output file"""

    @staticmethod
    def get_current():
        data = local()
        if not hasattr(data, 'model_logger'):
            data.model_logger = ModelLogger()
        return data.model_logger

    def __init__(self):
        self._lines = []

    def clear(self):
        self._lines = []

    def add(self, line):
        self._lines.append(line)

    def write(self, log_path):
        with open(log_path, 'w') as f:
            for line in self._lines:
                f.write(line + '\n')
