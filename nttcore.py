import numpy as np
from scipy.fft import fft, ifft, rfft, irfft

import itertools


def reciprocal(n, mod):
    x, y = mod, n
    a, b = 0, 1
    while y != 0:
        a, b = b, a - x // y * b
        x, y = y, x % y
    if x == 1:
        return a % mod
    else:
        print("Error-----------")


def unique_prime_factors(n):
    if n < 1:
        raise ValueError()
    result = []
    i = 2
    end = sqrt(n)
    while i <= end:
        if n % i == 0:
            n //= i
            result.append(i)
            while n % i == 0:
                n //= i
            end = sqrt(n)
        i += 1
    if n > 1:
        result.append(n)
    return result


def is_primitive_root(val, degree, mod):
    if not (0 <= val < mod):
        raise ValueError()
    if not (1 <= degree < mod):
        raise ValueError()
    pf = unique_prime_factors(degree)
    return pow(val, degree, mod) == 1 and all(
        (pow(val, degree // p, mod) != 1) for p in pf
    )


def is_generator(val, totient, mod):
    if not (0 <= val < mod):
        raise ValueError()
    if not (1 <= totient < mod):
        raise ValueError()
    pf = unique_prime_factors(totient)
    return pow(val, totient, mod) == 1 and all(
        (pow(val, totient // p, mod) != 1) for p in pf
    )


def sqrt(n):
    if n < 0:
        raise ValueError()
    i = 1
    while i * i <= n:
        i *= 2
    result = 0
    while i > 0:
        if (result + i) ** 2 <= n:
            result += i
        i //= 2
    return result


def is_prime(n):
    if n <= 1:
        raise ValueError()
    return all((n % i != 0) for i in range(2, sqrt(n) + 1))


def find_modulus(veclen, minimum):
    if veclen < 1 or minimum < 1:
        raise ValueError()
    start = (minimum - 1 + veclen - 1) // veclen
    for i in itertools.count(max(start, 1)):
        n = i * veclen + 1
        assert n >= minimum
        if is_prime(n):
            return n


def find_generator(totient, mod):
    if not (1 <= totient < mod):
        raise ValueError()
    for i in range(1, mod):
        if is_generator(i, totient, mod):
            return i
    raise ValueError("No generator exists")


def find_primitive_root(degree, totient, mod):
    if not (1 <= degree <= totient < mod):
        raise ValueError()
    if totient % degree != 0:
        raise ValueError()
    gen = find_generator(totient, mod)
    root = pow(gen, totient // degree, mod)
    assert 0 <= root < mod
    return root


def fftmul(x, y):
    L = len(x) + len(y) - 1
    A = fft(x, L)
    B = fft(y, L)
    C = A * B
    D = ifft(C)
    D = abs(D)
    # D = D.astype(int)
    return D


def DFTslow(x):
    """Compute the discrete Fourier Transform of the 1D array x"""
    x = np.asarray(x, dtype=float)
    N = x.shape[0]
    n = np.arange(N)
    k = n.reshape((N, 1))
    M = np.exp(-2j * np.pi * k * n / N)
    return np.dot(M, x)


def ntt(x, r, mod):
    x = np.asarray(x, dtype=int)
    N = x.shape[0]
    n = np.arange(N)
    k = n.reshape((N, 1))
    # M = np.exp(-2j * np.pi * k * n / N)
    # M=np.zeros([N,N])
    coef = n * k
    M = np.power(1, n * k)
    for i in range(0, N):
        for j in range(0, N):
            M[i, j] = pow(int(r), int(coef[i, j]), mod)

    # M=M%mod
    sm = np.dot(M, x)
    rowsums = sm.sum(axis=1)
    # quotient=np.mod(rowsums,mod)
    quotient = rowsums % mod
    return quotient.reshape(N, 1)


def intt(x, r, n, N, mod):
    nr = pow(r, N - 2, N)
    tm = ntt(x, nr, N)
    sc = pow(n, N - 2, N)
    return (sc * tm) % N


def circular_convolve(vec0, vec1, n, mod):
    maxval = max(val for val in itertools.chain(vec0, vec1))
    minmod = maxval**2 * len(vec0) + 1
    N = find_modulus(n, int(minmod))
    r = find_primitive_root(n, N - 1, N)
    temp0 = ntt(vec0, r, N)
    temp1 = ntt(vec1, r, N)
    temp2 = (temp0 * temp1) % N
    return intt(temp2, r, n, N, mod)


def nttL(x, r, mod, L):
    n = len(x)
    if L > n:
        d = L - n
        x = np.pad(x, (0, d), "constant")
        n = len(x)

    x = x.reshape((n, 1))
    x = np.asarray(x, dtype=int)
    row = np.arange(n)
    col = row.reshape((n, 1))
    coef = row * col
    M = np.ones((n, n), dtype=int)

    for i in range(0, n):
        for j in range(0, n):
            M[i, j] = pow(int(r), int(coef[i, j]), mod)

    tmp = np.dot(M, x)
    rowsums = tmp.sum(axis=1)
    print(rowsums)
    quotient = rowsums % mod
    print("n: " + str(n))
    return quotient.reshape(n, 1)
