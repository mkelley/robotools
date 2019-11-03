from contextlib import contextmanager

@contextmanager
def file_lock(filename):
    import os
    lock = open(filename, 'w')
    yield
    lock.close()
    os.unlink(filename)
