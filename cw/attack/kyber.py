import numpy as np
import kyber_ntt as ntt

n = 256
q = 3329

POLY_BYTES = 384
KYBER768_DU = 10
KYBER768_POLYCOMPRESSEDBYTES_DU = 320
KYBER768_ETA1 = 2

def central_reduce(x):
    t = x % q
    return t - (t > q//2)*q


def bytes_to_bits(input_bytes):
    bit_string = ''.join(format(byte, '08b')[::-1] for byte in input_bytes)
    return list(map(int, list(bit_string)))


def decode(byte_arr, l):
    poly = np.zeros(n, dtype='uint16')
    B = bytes_to_bits(byte_arr)
    for i in range(n):
        poly[i] = sum([B[i * l + j] * 2 ** j for j in range(l)])
    return poly


def decode_kyber768_du(byte_arr):
    return decode(byte_arr, KYBER768_DU)


def decode_sk(byte_arr):
    return decode(byte_arr, 12) #Kyber's q is 12-bit


def bitstring_to_bytes(s):
    return bytes([int(s[i:i + 8][::-1], 2) for i in range(0, len(s), 8)])


def encode(poly, l):
    bit_string = ''.join(format(c, f'0{l}b')[::-1] for c in poly)
    return bitstring_to_bytes(bit_string)


def compress(poly, d):
    return (poly * ((1 << d))).astype('uint16') % (1 << d)


def round_half_up(x):
    out = x.copy()
    mask = (out >= 0)
    np.add(out, 0.5, where=mask, out=out)
    np.floor(out, where=mask, out=out)
    np.invert(mask, out=mask)
    np.subtract(out, 0.5, where=mask, out=out)
    np.ceil(out, where=mask, out=out)
    return out


def decompress(poly, d):
    return round_half_up(poly * (q / (1 << d))).astype('uint16')


def kyber768_decompress_du(poly):
    return decompress(poly, KYBER768_DU)


def kyber768_poly_unpack_decompress_du(byte_arr):
    assert (len(byte_arr)) == KYBER768_POLYCOMPRESSEDBYTES_DU
    poly = decode_kyber768_du(byte_arr)
    poly_d = kyber768_decompress_du(poly)
    return ntt.ntt(poly_d)