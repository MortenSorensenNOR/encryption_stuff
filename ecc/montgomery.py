import math
import random
import numpy as np

cycles = 0

def montgom_mult(a, b, p, r_bits):
    """
    Compute a * b * R^(-1) mod p
    """
    # precompute p' where p * p' = -1 mod R
    p_prime = pow(-p, -1, 2**r_bits)

    # standard montgomery
    t = a * b
    m = ((t % (2**r_bits)) * p_prime) % (2**r_bits)
    t = (t + m * p) >> r_bits

    if t >= p:
        t -= p

    return t

def mul_word_by_full(word, full, word_bits=32, total_bits=256):
    global cycles

    """
    Calculate word x full-width using only word-by-word multiplies
    """
    n_words = total_bits // word_bits 
    word_mask = (1 << word_bits) - 1

    full_words = [(full >> (i * word_bits)) & word_mask for i in range(n_words)]

    result = 0
    for i, fw in enumerate(full_words):
        partial = word * fw
        result += partial << (i * word_bits)

        cycles += 1

    return result

def montgom_mult_serial(a, b, p, r_bits, word_bits=16):
    global cycles

    n_words = r_bits // word_bits
    word_mask = (1 << word_bits) - 1

    # precompute p' for single word
    p_prime = pow(-p, -1, 2**word_bits) & word_mask

    # split "a" into words
    a_words = [(a >> (i * word_bits)) & word_mask for i in range(n_words)]

    # main multiplication loop
    t = 0
    for i in range(n_words):
        # add a[i] * b to accumulator
        t = t + mul_word_by_full(a_words[i], b)

        # compute reduction factor for this iteration
        m = ((t & word_mask) * p_prime) & word_mask
        cycles += 1

        # add m * p and shift
        t = (t + mul_word_by_full(m, p)) >> word_bits

    if t >= p:
        t -= p
        cycles += 1

    return t

def mod_add(a, b, p):
    global cycles
    s = a + b
    cycles += 1

    if s >= p:
        s -= p
        cycles += 1

    return s

def mod_sub(a, b, p):
    global cycles
    cycles += 1
    if a >= b:
        return a - b
    return p - (b - a)

def mont_square(a, p, r_bits):
    """Square in Montgomery domain"""
    return montgom_mult_serial(a, a, p, r_bits)

def print_cycles():
    global cycles
    print(f"Total cycles: {cycles}")
