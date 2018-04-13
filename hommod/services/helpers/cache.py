import logging
import datetime

import redis
import pickle

from dogpile.cache.util import function_key_generator

from hommod.models.error import ServiceError


__all__ = ['cache_manager']


_log = logging.getLogger(__name__)


class CacheManager(object):
    def __init__(self, redis_hostname=None, redis_port=None,
                 redis_db=None, expiration_time=None, lock_timeout=None):
        self.redis_hostname = redis_hostname
        self.redis_port = redis_port
        self.redis_db = redis_db
        self.expiration_time = expiration_time
        self.lock_timeout = lock_timeout

        self._enabled = True

    def disable(self):
        self._enabled = False

    def enable(self):
        self._enabled = True

    def _get_redis(self):
        if self.redis_hostname is None:
            raise ServiceError("redis hostname is not set")

        if self.redis_port is None:
            raise ServiceError("redis port is not set")

        if self.redis_db is None:
            raise ServiceError("redis db is not set")

        return redis.StrictRedis(host=self.redis_hostname, port=self.redis_port, db=self.redis_db)

    def _get_key(self, f, args, kwargs):
        key = function_key_generator(None, f)(*args, **kwargs)
        return key

    def _get_lock_name(self, f, args, kwargs):
        return 'lock_%s' % self._get_key(f, args, kwargs)

    def _get_time_key(self, f, args, kwargs):
        return 'time_%s' % self._get_key(f, args, kwargs)

    def _get_value(self, r, f, args, kwargs):
        key = self._get_key(f, args, kwargs)

        time_str = r.get(self._get_time_key(f, args, kwargs))
        if time_str is None:
            return None

        time = datetime.datetime.strptime(time_str.decode('ascii'), "%Y-%m-%dT%H:%M:%S")
        now = datetime.datetime.now()
        delta = now - time
        _log.debug("delta for value is ({}) {}, thus {} seconds".format(type(delta), delta, delta.seconds))
        if delta.seconds > self.expiration_time:
            _log.debug("value has expired for {}.{}".format(f.__module__, f.__name__))
            return None

        value = r.get(key)
        if value is None:
            return None
        return pickle.loads(value)

    def _set_value(self, r, f, args, kwargs, value):
        key = self._get_key(f, args, kwargs)

        now = datetime.datetime.now()
        r.set(self._get_time_key(f, args, kwargs), now.strftime("%Y-%m-%dT%H:%M:%S").encode('ascii'))

        r.set(key, pickle.dumps(value))

    def delete(self, f, *args, **kwargs):
        r = self._get_redis()
        key = self._get_key(f, args, kwargs)

        r.delete(key)

    def cache(self):
        def wrapped(f):
            def new_f(*args, **kwargs):

                if not self._enabled:
                    return f(*args, **kwargs)

                _log.debug("calling new_f")

                r = self._get_redis()
                lock = redis.lock.Lock(r, self._get_lock_name(f, args, kwargs),
                                       blocking_timeout=self.lock_timeout)

                if lock.acquire():
                    _log.debug('lock success for {}.{}'.format(f.__module__, f.__name__))
                    try:
                        value = self._get_value(r, f, args, kwargs)
                        if value is not None:

                            _log.debug('returning old value for {}.{}'.format(f.__module__, f.__name__))
                            return value
                        else:
                            _log.debug('setting new value for {}.{}'.format(f.__module__, f.__name__))

                            value = f(*args, **kwargs)
                            self._set_value(r, f, args, kwargs, value)
                            return value
                    finally:
                        lock.release()
                else:
                    _log.debug('lock failed for {}.{}'.format(f.__module__, f.__name__))

                    value = self._get_value(r, f, args, kwargs)
                    if value is not None:
                        _log.debug('returning old value for {}.{}'.format(f.__module__, f.__name__))
                        return value
                    else:
                        _log.debug('computing value for {}.{}'.format(f.__module__, f.__name__))
                        return f(*args, **kwargs)

            new_f.__name__ = f.__name__
            new_f.__module__ = f.__module__
            return new_f
        return wrapped


cache_manager = CacheManager()
