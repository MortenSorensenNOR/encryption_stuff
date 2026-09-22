import math

MODULUS = 17
GEN     = 9

def brv(x, n):
    return int(''.join(reversed(bin(x)[2:].zfill(n))), 2)

def ntt_iter(a, gen=GEN, modulus=MODULUS):
    deg_d = len(a)

    stride = 1
    nbits = int(math.log2(deg_d))
    res   = [a[brv(i, nbits)] for i in range(deg_d)]

    # precompute generators used in different stages of recursion
    gens = [pow(gen, pow(2, i), modulus) for i in range(nbits)]
    # first layer uses the lowest root of unity, i.e. last one
    gen_ptr = len(gens) - 1

    while stride < deg_d:
        for start in range(0, deg_d, stride*2):
            for i in range(start, start + stride):
                # compute the omega multiplier - j = i - start
                zp = pow(gens[gen_ptr], i - start, modulus)

                # Cooley-Tukey butterfly
                a = res[i]
                b = res[i+stride]
                res[i]        = (a + zp * b) % modulus
                res[i+stride] = (a - zp * b) % modulus

        # grow stride
        stride <<= 1
        # move to the next root of unity
        gen_ptr -= 1
    
    return res

def intt_iter(a, gen=GEN, modulus=MODULUS):
    deg_d = len(a)

    stride = deg_d // 2

    # shuffle input array in bit-reversal order
    nbits  = int(math.log2(deg_d))
    res = a[:]

    # precompute inverse generators used in different stages of recursion
    gen = pow(gen, -1, modulus)
    gens = [pow(gen, pow(2, i), modulus) for i in range(nbits)]
    gen_ptr = 0

    while stride > 0:
        for start in range(0, deg_d, stride*2):
            for i in range(start, start + stride):
                zp = pow(gens[gen_ptr], i - start, modulus)

                # Gentleman-Sande butterfly
                a = res[i]
                b = res[i+stride]
                res[i] = (a + b) % modulus
                res[i+stride] = ((a - b) * zp) % modulus

        stride >>= 1
        gen_ptr += 1

    scaler = pow(deg_d, -1, modulus)
    return [(res[brv(i, nbits)] * scaler) % modulus for i in range(deg_d)]

a     = [1, 2, 3, 4, 0, 0, 0, 0]
b     = [4, 2, 1, 8, 0, 0, 0, 0]
a_hat = ntt_iter(a)
b_hat = ntt_iter(b)
c_hat = [a_hat[i] * b_hat[i] % MODULUS for i in range(len(a))] 

c = intt_iter(c_hat)
print(a, b, c)
