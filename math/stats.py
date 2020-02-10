"""Helpful statistical functions."""

from statistics import StatisticsError, mean
from math import comb, exp, sqrt, pi, erf
from functools import partial

def binomial_prob(p, n, k):
    """Finds the probability that the event with probability 'p' will
    occur a given number of times in a given number of trials."""
    q = 1 - p
    return comb(n, k) * p**k * q**(n - k)

def quantile(n):
    """Returns a function that calculates the n-quantiles of a data set."""
    def inner(data, k):
        if not 0 < k < n:
            raise StatisticsError
        data = sorted(data)
        rank = len(data) * (k/n)
        if float.is_integer(rank):
            return mean((data[int(rank) - 1], data[int(rank)]))
        else:
            return data[int(rank)]
    return inner

quartile = quantile(4)
quartile.__doc__ = "Returns the `k`th quartile of a data set."

percentile = quantile(100)
percentile.__doc__ = "Returns the `k`th percentile of a data set."""

class Distribution:
    """Base class for probability distributions."""

    def __init__(self, mu, sigma):
        """Initialize self. See help(type(self)) for accurate signature."""
        self.mu = mu
        self.sigma = sigma

    def pdf(self, x):
        """Probability density function."""
        pass

    def cdf(self, x):
        """Cumulative distribution function."""
        pass

class NormalDistribution(Distribution):
    """The normal distribution."""

    def pdf(self, x):
        """Probability density function."""
        return exp((-x - self.mu)**2 / (2 * self.sigma**2)) / \
               sqrt(2*pi * self.sigma**2)

    def cdf(self, x):
        """Cumulative distribution function."""
        return (1 + erf((x - self.mu) / (self.sigma * sqrt(2)))) / 2

std_normal = NormalDistribution(0, 1)

class StudentDistribution(Distribution):
    """Student's t-distribution."""

    def __init__(self, mu, sigma, v):
        """Initialize self. See help(type(self)) for accurate signature."""
        super().__init__(mu, sigma)
        self.v = v

    def pdf(self, x):
        """Probability density function."""
        pass

    def cdf(self, x):
        """Cumulative distribution function."""
        pass
