"""Various tools for the representation and study of sequences."""

from functools import partial, reduce
from operator import ne, add
from math import gcd
from fractions import Fraction

from sympy.ntheory import totient

from functions import lcm, extended_gcd

def extended_fibonacci(x, y):
    """An extension of the Fibonacci sequence to any two integers."""
    yield x
    yield y
    while True:
        z = x + y
        yield z
        x, y = y, z

fibonacci = partial(extended_fibonacci, 0, 1)
fibonacci.__doc__ = "The full (infinite) Fibonacci sequence."

lucas = partial(extended_fibonacci, 2, 1)
lucas.__doc__ = "The full (infinite) Lucas sequence."

def collatz(n):
    """The Collatz sequence beginning with 'n'."""
    while n != 1:
        yield n
        n = n // 2 if n % 2 == 0 else 3*n + 1
    yield 1

def shortcut_collatz(n):
    """The Collatz sequence beginning with 'n' from the shortcut
    Collatz map."""
    while n > 1:
        yield n
        if n % 2:
            n = 3*n + 1
        n //= 2
    yield 1

def reductive_lcm(n):
    """A sequence such that the 'n'th term of the sequence is the least
    common multiple of all positive integers less than or equal to 'n'."""
    for i in range(2, n+1):
        yield lcm(range(1, i+1))

def least_non_divisors(a, b):
    """Yields the least integers 'm' such that n % m != 0, for all n
    in [a, b)."""
    def inner(n):
        m = 1
        while not n % m:
            m += 1
        return m
    for n in range(a, b):
        yield inner(n)

def hamming_weight(n, subseq=None):
    """Returns the Hamming weights of all non-negative integers in the
    interval [0, 2**n).

    The terms of this sequence are generated via a recurrence relation.
    If 'subseq' is passed all the terms in [0, 2**(n-1)), the rest of
    of the terms can be generated more efficiently.
    """
    weights = subseq if subseq else [0]
    while len(weights) < 2**n:
        for i in weights.copy():
            weights.append(i + 1)
    return weights

def bit_inversions(n):
    """Yields the number of inversions required to shift all the ones
    in a binary number to the right side, for all non-negative integers
    less than 2**n."""
    yield 0
    terms = [0]
    weights = [0]
    for i in range(n):
        weights = hamming_weight(i, weights)
        for j in range(len(weights)):
            term = terms[j] + weights[-(j + 1)]
            yield term
            terms.append(term)

def farey_length(n):
    """Yields the successive lengths of Farey sequences."""
    length = 1
    for i in range(n):
        yield length
        length += totient(i + 1)

def farey(n, first=(0, 1), last=(1, 1)):
    """Yields the terms of the nth Farey sequence from 'first' to
    'last', inclusive.

    This function utilizes the algorithm described in the Wikipedia
    article on Farey sequences to generate every term except the first
    and second.
    """
    # reduces 'first' and 'last'
    r, s = Fraction(*first), Fraction(*last)
    a, b, last = r.numerator, r.denominator, (s.numerator, s.denominator)
    # finds second term
    x, y = extended_gcd(b, -a)[1:]
    r = (n - y) // b
    c, d = r*a + x, r*b + y
    while (a, b) != last:
        yield (a, b)
        r = (n + b) // d
        a, b, c, d = c, d, r*c - a, r*d - b
    yield (a, b)

def perfect_squares(a, b):
    """Yields the perfect squares in the interval [a**2, b**2)."""
    square = a**2
    while square < b**2:
        yield square
        a += 1
        square += 2*a - 1

def gcd_sequence(n):
    """The sequence of gcd(n, k) for all 'k' less than 'n'."""
    for k in range(1, n+1):
        yield gcd(n, k)

def totatives(n):
    """Yields the integers 'a' such that a <= n and gcd(a, n) = 1."""
    for k in range(1, n+1):
        if gcd(n, k) == 1:
            yield k
