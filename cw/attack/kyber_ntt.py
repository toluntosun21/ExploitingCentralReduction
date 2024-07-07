import numpy as np
import math

########################### KYBER SIMPLE NTT #############################


q = 3329
root = 17
root_inv = pow(root, -1, q)

def reverse_bits(n, width=7):    
    b = '{:0{width}b}'.format(n, width=width)
    return int(b[::-1], 2)


__n = 128
n = 256

bit_reverse_table = [reverse_bits(i) for i in range(__n)]


__ntt_mat = np.zeros((__n,__n),dtype='int16')
for i in range(__n):
    for j in range(__n):
        __ntt_mat[i,j] = pow(root, (2*bit_reverse_table[i] + 1)*j, q)

__ntt_mat_inv = np.zeros((__n,__n), dtype='int16')
inv2 = pow(__n, -1, q)
for i in range(__n):
    for j in range(__n):
        __ntt_mat_inv[j,i] = (pow(root_inv, (2*bit_reverse_table[i] + 1)*j, q) * inv2) % q


def ntt(a, central_red=True):
    temp = np.empty((256), dtype='int16')
    temp[::2] = (np.matmul(__ntt_mat.astype('int64'), a[::2].astype('int64')) % q).astype('int16')
    temp[1::2] = (np.matmul(__ntt_mat.astype('int64'), a[1::2].astype('int64')) % q).astype('int16')
    if central_red:
        return temp - (temp > q//2)*q # centralize
    else:
        return temp


def ntt_inv(a, central_red=True):
    temp = np.empty((256), dtype='int16')
    temp[::2] = (np.matmul(__ntt_mat_inv.astype('int64'), a[::2].astype('int64')) % q).astype('int16')
    temp[1::2] = (np.matmul(__ntt_mat_inv.astype('int64'), a[1::2].astype('int64')) % q).astype('int16')
    if central_red:
        return temp - (temp > q//2)*q # centralize
    else:
        return temp


__ntt_mat = np.zeros((__n,__n),dtype='uint16')
for i in range(__n):
    for j in range(__n):
        __ntt_mat[i,j] = pow(root, (2*bit_reverse_table[i] + 1)*j, q)

__ntt_mat_inv = np.zeros((__n,__n), dtype='uint16')
inv2 = pow(__n, -1, q)
for i in range(__n):
    for j in range(__n):
        __ntt_mat_inv[j,i] = (pow(root_inv, (2*bit_reverse_table[i] + 1)*j, q) * inv2) % q


def ntt(a, central_red=False):
    temp = np.empty((256), dtype='uint16')
    temp[::2] = (np.matmul(__ntt_mat.astype('uint64'), a[::2].astype('int64')) % q).astype('uint16')
    temp[1::2] = (np.matmul(__ntt_mat.astype('uint64'), a[1::2].astype('int64')) % q).astype('uint16')
    if central_red:
        return (temp - (temp > q//2)*q).astype('int16') # centralize
    else:
        return temp


def ntt_inv(a, central_red=False):
    temp = np.empty((256), dtype='uint16')
    temp[::2] = (np.matmul(__ntt_mat_inv.astype('uint64'), a[::2].astype('int64')) % q).astype('uint16')
    temp[1::2] = (np.matmul(__ntt_mat_inv.astype('uint64'), a[1::2].astype('int64')) % q).astype('uint16')
    if central_red:
        return (temp - (temp > q//2)*q).astype('int16') # centralize
    else:
        return temp



####################################### PLANTARD #####################################################



# Plantard arithmetic models the implementation in https://eprint.iacr.org/2022/112.pdf
# with improved plantard reduction implementation https://eprint.iacr.org/2022/956.pdf
# included in pqm4 library, commit 3743a66

q = 3329
qa = 26632
qinv = 0x6ba8f301

inv2 = pow(2,-1,q)


zetas_plant = [21932846, 3562152210, 752167598, 3417653460, 2112004045, 932791035, 2951903026, 1419184148, 1817845876, 3434425636, 4233039261, 300609006, 975366560, 2781600929, 3889854731, 3935010590, 2197155094, 2130066389, 3598276897, 2308109491, 2382939200, 1228239371, 1884934581, 3466679822, 1211467195, 2977706375, 3144137970, 3080919767, 945692709, 3015121229, 345764865, 826997308, 2043625172, 2964804700, 2628071007, 4154339049, 483812778, 3288636719, 2696449880, 2122325384, 1371447954, 411563403, 3577634219, 976656727, 2708061387, 723783916, 3181552825, 3346694253, 3617629408, 1408862808, 519937465, 1323711759, 1474661346, 2773859924, 3580214553, 1143088323, 2221668274, 1563682897, 2417773720, 1327582262, 2722253228, 3786641338, 1141798155, 2779020594]


_zetas_temp = []
for zeta in zetas_plant:
    _zetas_temp.append(zeta)
    _zetas_temp.append(0xFFFFFFFF - zeta)
zetas_plant_np = np.array(_zetas_temp, dtype='uint32').astype('int32')



def plant_red(in_, factor=qinv):
    a1 = in_.astype('int32')

    t0 = (((np.uint32(factor).astype('int32').astype('int64') * a1.astype('int64')).astype('uint64') >> np.uint64(16))).astype('uint16').astype('int16').astype('int32')
    t1 = ((np.int32(qa) + t0 * np.int32(q)).astype('uint32') >> 16).astype('uint16')
    return t1.astype('int16')



def plant_red_zetas(in_, zetas=zetas_plant_np):
    a1 = in_.astype('int32')

    t0 = (((zetas.astype('int32').astype('int64') * a1.astype('int64')).astype('uint64') >> np.uint64(16))).astype('uint16').astype('int16').astype('int32')
    t1 = ((np.int32(qa) + t0 * np.int32(q)).astype('uint32') >> 16).astype('uint16')
    return t1.astype('int16')


#poly_frombytes_mul_32_16
def basemul_plant_low(a, b, frame=range(0,__n)):
    t0 = plant_red_zetas(b[...,1::2].astype(np.int16), zetas_plant_np[frame])
    t = a[...,::2].astype(np.int16).astype(np.int32) * b[...,::2].astype(np.int16).astype(np.int32) + a[...,1::2].astype(np.int16).astype(np.int32) * t0.astype(np.int32)
    return plant_red(t)


#poly_frombytes_mul_32_16
def basemul_plant_high(a, b):
    t = a[...,::2].astype(np.int16).astype(np.int32) * b[...,1::2].astype(np.int16).astype(np.int32) + a[...,1::2].astype(np.int16).astype(np.int32) * b[...,::2].astype(np.int16).astype(np.int32)
    return plant_red(t)


#poly_frombytes_mul_32_16
def basemul_plant(a, b, frame=range(0,__n)):
    r = np.empty((a.shape[0]), dtype='int16')
    t0 = plant_red_zetas(b[1::2].astype(np.int16), zetas_plant_np[frame])
    t = a[::2].astype(np.int16).astype(np.int32) * b[::2].astype(np.int16).astype(np.int32) + a[1::2].astype(np.int16).astype(np.int32) * t0.astype(np.int32)
    r[::2] = plant_red(t)
    t = a[::2].astype(np.int16).astype(np.int32) * b[1::2].astype(np.int16).astype(np.int32) + a[1::2].astype(np.int16).astype(np.int32) * b[::2].astype(np.int16).astype(np.int32)
    r[1::2] = plant_red(t)
    return r


#poly_frombytes_mul_32_16
def basemul_plant(a, b, frame=range(0,__n)):
    r = np.empty(shape=a.shape, dtype='int16')
    r[...,::2] = basemul_plant_low(a, b, frame=frame)
    r[...,1::2] = basemul_plant_high(a, b)
    return r



####################################### MONTGOMERY #####################################################

# based on the implementation in https://github.com/uclcrypto/pqm4_masked
# kyber768
# commit 5fe90ba

mu = 3327

zetas_monty = [2226, 430,  555,  843,  2078, 871,  1550, 105,  422,  587,  177, 3094, 3038, 2869, 1574, 1653, 3083, 778,  1159, 3182, 2552, 1483, 2727, 1119, 1739, 644,  2457, 349,  418,  329,  3173, 3254, 817, 1097, 603,  610,  1322, 2044, 1864, 384,  2114, 3193, 1218, 1994, 2455, 220,  2142, 1670, 2144, 1799, 2051, 794,  1819, 2475, 2459, 478,  3221, 3021, 996,  991,  958,  1869, 1522, 1628]


__zetas_monty_temp = []
for zeta in zetas_monty:
    __zetas_monty_temp.append(zeta)
    __zetas_monty_temp.append(-zeta)
zetas_monty_np = np.array(__zetas_monty_temp, dtype='int16')



def montgomery(in_):
    t = in_.astype(np.int32)
    t1 = t.astype(np.int16).astype(np.int32) * np.array([mu], dtype=np.int32)
    t2 = t1.astype(np.int16).astype(np.int32) * np.array([q], dtype=np.int16).astype(np.int32)
    t3 = ((t + t2).astype(np.int32))
    t4 = (t3 >> 16) & 0xFFFF
    return t4.astype(np.int16)



def basemul_monty_low(a, b, frame=range(0,__n)):
    t0 = montgomery(a[...,1::2].astype(np.int16).astype(np.int32) * b[...,1::2].astype(np.int16).astype(np.int32))
    t = a[...,::2].astype(np.int16).astype(np.int32) * b[...,::2].astype(np.int16).astype(np.int32) + zetas_monty_np[frame] * t0.astype(np.int16).astype(np.int32)
    return montgomery(t)



def basemul_monty_high(a, b):
    t = a[...,::2].astype(np.int16).astype(np.int32) * b[...,1::2].astype(np.int16).astype(np.int32) + a[...,1::2].astype(np.int16).astype(np.int32) * b[...,::2].astype(np.int16).astype(np.int32)
    return montgomery(t)



def basemul_monty(a, b, frame=range(0,__n)):
    r = np.empty(shape=a.shape, dtype='int16')
    r[...,::2] = basemul_monty_low(a, b, frame=frame)
    r[...,1::2] = basemul_monty_high(a, b)
    return r