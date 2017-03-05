import time


def timeit(func):
    def wrapper(*args, **kwargs):
        start = time.time()
        func(*args, **kwargs)
        total = time.time() - start
        print("total ex time: {}".format(total))
    return wrapper
