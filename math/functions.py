"""Various mathematical functions not found in the standard library."""

from collections import namedtuple
from collections.abc import Iterable
from math import gcd
from functools import reduce, partial, lru_cache
from operator import mul
from fractions import Fraction

import sympy.ntheory as nthry

def lcm(*args):
    """Finds the least common multiple of an arbitrary number of
    integers.

    Always returns a positive integer, even when some or all of the
    arguments are negative.
    """
    if isinstance(args[0], int):
        iterable = args
    elif isinstance(args[0], Iterable):
        iterable = args[0]
    return reduce(lambda a, b: (a // gcd(a, b)) * b, iterable, 1)

def double_factorial(n):
    """Returns the double factorial of 'n'."""
    if n in range(-1, 2):
        return 1
    else:
        return reduce(mul, (k for k in range(n, 1, -2)))

def multifactorial(n, k):
    """A generalization of the factorial function."""
    pass

def catalan(n):
    """Returns the 'n'th Catalan number."""
    pass

def primorial_totient(n):
    """Returns the totient of a given primorial."""
    return reduce(mul, (p-1 for p in nthry.primerange(2, nthry.prime(n + 1))))

@lru_cache(None)
def fusc(n):
    """Returns a given term of Stern's diatomic sequence."""
    if n == 0:
        return 0
    elif n & (n - 1) == 0:
        return 1
    else:
        if not n % 2:
            return fusc(n // 2)
        else:
            return fusc(n // 2) + fusc(n // 2 + 1)

def convolution(f, g):
    """Returns the Dirichlet convolution of two functions."""
    def func(n):
        return sum(f(d)*g(n//d) for d in nthry.divisors(n, True))
    return func

divisor_sum = partial(convolution, g=lambda x: 1)
divisor_sum.__doc__ = "The sum of a function applied to all of a number's " \
                      "divisors."

mobius_inverse = partial(convolution, g=nthry.mobius)
mobius_inverse.__doc__ = "The Mobius inverse of a function."

gcd_composition = partial(convolution, g=nthry.totient)
gcd_composition.__doc__ = "For a function 'f', returns a function 'h' such " \
                          "that 'h(n) == sum(f(gcd(n, k)) for k in " \
                          "range(1, n+1).'"

# Handle cases where f(1) should equal zero but is a value extremely
# close to zero instead.
def dirichlet_inverse(f):
    if not f(1):
        raise ValueError("Function 'f' is not invertible, as f(1) == 0.")
    @lru_cache(None)
    def g(n):
        divs = nthry.divisors(n)[:-1]
        x = -sum(f(n//d)*g(d) for d in divs) if divs else 1
        return x / f(1)
    return g

def _upsilon(n):
    """Yields the numbers unrelated to n."""
    for k in range(1, n):
        if 1 < gcd(n, k) < k:
            yield k

def unrelated(n):
    """The amount of numbers less than n that are neither coprime to
    nor divisors of n."""
    return n - nthry.totient(n) - nthry.divisor_count(n) + 1

def digamma(n, a, b):
    """A generalization of sum of prime divisors."""
    factors = nthry.factorint(n)
    return sum(p**b for p in factors if factors[p] > a)

def cfract(a, b):
    """Returns the continued fraction representation of a/b."""
    return list(nthry.continued_fraction_iterator(Fraction(a, b)))

Euclidean = namedtuple("Euclidean", ["gcd", "x", "y"])

def extended_gcd(a, b):
    """Implementation of the extended Euclidean algorithm.

    For ax + by = gcd(a, b) = d, returns the 'namedtuple' (d, x, y),
    accessible via the fields "gcd", "x", and "y", respectively.
    """
    x1, y1, x2, y2 = 1, 0, 0, 1
    while b:
        div, mod = divmod(a, b)
        a, b = b, mod
        x1, y1 = y1, x1 - div*y1
        x2, y2 = y2, x2 - div*y2
    sign = 1 if a > 0 else -1
    return Euclidean(sign * a, sign * x1, sign * x2)
