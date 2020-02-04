"""Contains tools that approximate various functions involving limits."""

from math import log, factorial

from pyutils.mathematics import EULER_MASC

def derivative(func, x, delta_x):
    """Approximates the derivative of a function at 'x' using its limit
    definition."""
    dividend = func(x + delta_x) - func(x)
    deriv = dividend / delta_x
    return deriv

def def_integral(func, a, b, n):
    """Approximates the definite integral of a function over the interval
    [a, b] via Riemann sums."""
    delta_x = (b - a) / n
    c = a
    total = 0
    for i in range(n + 1):
        total += func(c) * delta_x
        c += delta_x
    return total

# Figure out how to pass in initial terms. (Perhaps lambda closures?)
def taylor_series(func, start):
    """Generates a function that finds the partial sum of a Taylor series
    with a given number of terms.

    'func' should have three parameters: the first being the center of
    the series, the second being the main argument, and the third being
    the number of terms to iterate over (the variable at the top of the
    summation symbol). 'start' is the initial index (the variable at the
    bottom of the summation symbol).
    """
    def partial_sum(c, x, n):
        """Finds the partial sum of 'n' terms of this Taylor series,
        which is centered at 'c' and evaluated at 'x'."""
        total = 0
        for i in range(start, start + n + 1):
            total += func(c, x, i)
        return total
    return partial_sum

def li(x, n):
    """Approximates the logarithmic integral of 'x' using an infinite
    series approximation discovered by Ramanujan."""
    u = log(x)
    result = EULER_MASC + log(u)
    for i in range(1, n+1):
        result += u**i / (i * factorial(i))
    return result

    

    
