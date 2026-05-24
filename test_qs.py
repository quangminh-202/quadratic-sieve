import math

from qs_course.quadratic_sieve import factor_integer, quadratic_sieve


def test_known_small_semiprime():
    n = 103 * 107
    g = quadratic_sieve(n, B=100, verbose=False)
    assert 1 < g < n
    assert n % g == 0


def test_recursive_factorization():
    n = 1009 * 1013
    factors = factor_integer(n, verbose=False, B=300)
    assert math.prod(factors) == n
    assert factors == [1009, 1013]

