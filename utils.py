from functools import wraps
from itertools import starmap, repeat
from string import digits
from random import choices

def memoized(f):
    memo = {}
    @wraps(f)
    def _(*args):
        if args not in memo:
            memo[args] = f(*args)
        return memo[args]
    return _

def ensure_unique(iterable, max_skips=float('inf')):
    """Yields the unique elements of an iterable, skipping those that
    have been encountered before. If more than `max_skips` values have
    have skipped, returns."""
    old = set()
    iterator = iter(iterable)
    try:
        while True:
            i = 0
            while i <= max_skips:
                x = next(iterator)
                if x not in old:
                    yield x
                    old.add(x)
                    break
                i += 1
            else:
                break
    except StopIteration:
        pass

# from itertools recipes
def repeatfunc(func, times=None, *args):
    """Repeat calls to func with specified arguments.

    Example:  repeatfunc(random.random)
    """
    if times is None:
        return starmap(func, repeat(args))
    return starmap(func, repeat(args, times))

def random_digits(n):
    return ''.join(choices(digits, k=n))
