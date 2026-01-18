from montgomery import *

# Curve25519 parameters
# Equation: y² = x³ + 486662x² + x (Montgomery curve, B=1)
p = 2**255 - 19
A = 486662
A24 = (A + 2) // 4  # = 121666, used in the ladder step

# For Montgomery domain
r_bits = 256  # We use 256 bits even though p is 255 bits
R = 2**256
R_inv = pow(R, -1, p)
R_mod_p = R % p

def ladder_step_mont_domain(X1, Z1, X2, Z2, X_diff, A24_mont, p, r_bits):
    """
    Montgomery ladder step in Montgomery domain arithmetic

    input:
        R0 = (X1, Z1) - curve point
        R1 = (X2, Z2) - curve point + base
        X_diff        - x-coordinate of base point (R1 - R0)
        A24_mont      - curve constant (A + 2) / 4 in montgomery domain

    output:
        (X1_new, Z1_new) = 2 * R0
        (X2_new, Z2_new) = R0 + R1
    """

    # precomputations for addition
    U = mod_sub(X1, Z1, p)
    V = mod_add(X2, Z2, p)
    U = montgom_mult_serial(U, V, p, r_bits) # U = (X1 - Z1) * (X2 + Z2)

    V = mod_add(X1, Z1, p)
    W = mod_sub(X2, Z2, p)
    V = montgom_mult_serial(V, W, p, r_bits) # V = (X1 + Z1) * (X2 - Z2)

    # addition result
    add_sum = mod_add(U, V, p)
    add_dif = mod_sub(U, V, p)

    X2_new = mont_square(add_sum, p, r_bits)                # (U + V)^2
    Z2_new = mont_square(add_dif, p, r_bits)                # (U - V)^2
    Z2_new = montgom_mult_serial(X_diff, Z2_new, p, r_bits) # X_diff * (U - V)^2

    # precomputes for doubling
    S = mod_add(X1, Z1, p)
    S = mont_square(S, p, r_bits)   # D = (X1 + Z1)^2

    D = mod_sub(X1, Z1, p)
    D = mont_square(D, p, r_bits)   # D = (X1 - Z1)^2
    
    X1_new = montgom_mult_serial(S, D, p, r_bits) # X_dbl = S * D

    E = mod_sub(S, D, p) # E = S - D

    T = montgom_mult_serial(A24_mont, E, p, r_bits) # A24 * E
    T = mod_add(D, T, p)                            # D + A24 * E
    Z1_new = montgom_mult_serial(E, T, p, r_bits)   # Z_dbl = E * (D + A24 * E)

    return X1_new, Z1_new, X2_new, Z2_new

def field_inv_mont(z, p, r_bits):
    """
    Compute z^(-1) in Montgomery domain using Fermat's little theorem.
    z^(-1) = z^(p-2) mod p for p = 2^255 - 19
    
    Cost: 254 squarings + 11 multiplications
    """
    def M(x, y):
        return montgom_mult_serial(x, y, p, r_bits)
    
    def S(x):
        return M(x, x)
    
    def pow2(x, n):
        """Square n times: x^(2^n)"""
        for _ in range(n):
            x = S(x)
        return x
    
    # Build up powers using addition chain
    t0 = S(z)           # z^2
    t1 = S(t0)          # z^4
    t2 = S(t1)          # z^8
    t3 = M(t2, z)       # z^9
    t4 = M(t3, t0)      # z^11
    t5 = M(S(t4), t3)   # z^31 = z^(2^5 - 1)
    
    t6 = M(pow2(t5, 5), t5)      # z^(2^10 - 1)
    t7 = M(pow2(t6, 10), t6)     # z^(2^20 - 1)
    t8 = M(pow2(t7, 20), t7)     # z^(2^40 - 1)
    t9 = M(pow2(t8, 10), t6)     # z^(2^50 - 1)
    t10 = M(pow2(t9, 50), t9)    # z^(2^100 - 1)
    t11 = M(pow2(t10, 100), t10) # z^(2^200 - 1)
    t12 = M(pow2(t11, 50), t9)   # z^(2^250 - 1)
    
    return M(pow2(t12, 5), t4)   # z^(2^255 - 21)

def x25519_scalar_mult(k, X_base, p, r_bits):
    """
    Compute k * P given only x-coordinate of P

    input: 
        k      - scalar
        X_base - x-coordinate of base point

    output:
        x-coordinate of k * P
    """

    # should be precomputed
    R        = 2**r_bits
    A24      = 121666
    A24_mont = (A24 * R) % p

    # convert base point to Montogmery domain
    X_diff = (X_base * R) % p # constant throughout compute

    # initialization
    one_mont = R % p
    zero     = 0

    X1 = one_mont # R0.X = 1 (at infinity)
    Z1 = zero     # R0.Z = 0
    X2 = X_diff   # R1.X = X_base (in mont domain)
    Z2 = one_mont # R1.Z = 1

    # process bits from top to bottom
    for i in range(254, -1, -1):
        ki = (k >> i) & 1

        if ki == 1:
            X1, X2 = X2, X1
            Z1, Z2 = Z2, Z1

        X1, Z1, X2, Z2 = ladder_step_mont_domain(X1, Z1, X2, Z2, X_diff, A24_mont, p, r_bits)

        if ki == 1:
            X1, X2, = X2, X1
            Z1, Z2 = Z2, Z1

    # compute X1/Z1 and convert back from montgomery domain
    R_inv = pow(R, -1, p)
    Z1_inv_mont = field_inv_mont(Z1, p, r_bits)
    result_mont = montgom_mult_serial(X1, Z1_inv_mont, p, r_bits)
    result = montgom_mult_serial(result_mont, 1, p, r_bits)

    return result

# not sure what the functions really do for testing, but the input test vectors are from some 
# spec with a certain encoding, so i guess abide by that
def bytes_to_int_le(hex_str):
    """Convert little-endian hex string to integer"""
    byte_arr = bytes.fromhex(hex_str)
    return int.from_bytes(byte_arr, 'little')

def int_to_bytes_le(n, length=32):
    """Convert integer to little-endian hex string"""
    return n.to_bytes(length, 'little').hex()

def clamp_scalar(k):
    """
    X25519 scalar clamping as per RFC 7748
    - Clear bottom 3 bits
    - Clear top bit
    - Set second-to-top bit
    """
    k &= ~7  # Clear bits 0, 1, 2
    k &= ~(1 << 255)  # Clear bit 255
    k |= (1 << 254)   # Set bit 254
    return k

def test_x25519_rfc():
    global cycles

    p = 2**255 - 19
    r_bits = 256
    
    # RFC 7748 test vector - decode as little-endian
    k_hex = "a546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4"
    u_hex = "e6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c"
    expected_hex = "c3da55379de9c6908e94ea4df28d084f32eccf03491c71f754b4075577a28552"
    
    k = bytes_to_int_le(k_hex)
    u = bytes_to_int_le(u_hex)
    expected = bytes_to_int_le(expected_hex)
    
    # Apply scalar clamping
    k = clamp_scalar(k)
    
    # Also need to mask u-coordinate (clear top bit)
    u &= (1 << 255) - 1
    
    print(f"Scalar k: {hex(k)}")
    print(f"Input u:  {hex(u)}")
    print(f"Expected: {hex(expected)}")
    
    result = x25519_scalar_mult(k, u, p, r_bits)
    
    print(f"Result:   {hex(result)}")
    print(f"Match: {result == expected}")

    print_cycles()

test_x25519_rfc()
