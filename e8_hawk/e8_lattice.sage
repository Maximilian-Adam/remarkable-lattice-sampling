import os
import sys

try:
    E8_HAWK_DIR
except NameError:
    _candidates = [
        os.path.dirname(os.path.abspath(__file__)),
        os.path.dirname(os.path.abspath(sys.argv[-1])) if sys.argv else "",
        os.path.join(os.getcwd(), "e8_hawk"),
        os.getcwd(),
    ]
    E8_HAWK_DIR = next(d for d in _candidates if os.path.exists(os.path.join(d, "params.sage")))

load(os.path.join(E8_HAWK_DIR, "params.sage"))

from sage.all import QQ, ZZ, codes, matrix, vector
import itertools


# Construction-A model used throughout this project:
#   E8_CA = RM(1,3) + 2 Z^8.
#
# This is the parity model that makes the b=2 cosets explicit. A vector is in
# this E8 model iff it is integral and its coordinate parities form a codeword
# in RM(1,3). Half shifts appear when sampling h + 2E8 by the identity
#   D_{h + 2E8, s} = 2 * D_{E8 + h/2, s/2}.

RM13_CODE = codes.BinaryReedMullerCode(1, 3)
RM13_CODEWORDS = [tuple(ZZ(x) for x in cw) for cw in RM13_CODE]
RM13_CODEWORD_SET = set(tuple(int(x) for x in cw) for cw in RM13_CODEWORDS)

_I2 = [vector(ZZ, [2 if i == j else 0 for j in range(E8_BLOCK_DIM)]) for i in range(E8_BLOCK_DIM)]
_G = RM13_CODE.generator_matrix()
_GROWS = [vector(ZZ, [ZZ(x) for x in _G.row(i)]) for i in range(_G.nrows())]
E8_BASIS = matrix(ZZ, _I2 + _GROWS).row_module().basis_matrix()
E8_GRAM = E8_BASIS * E8_BASIS.transpose()


def vector_over(ring, v):
    return vector(ring, list(v))


def is_integral_vector(v):
    return all(QQ(x) in ZZ for x in v)


def parity_tuple(v):
    if not is_integral_vector(v):
        raise ValueError("parity is only defined for integral vectors")
    return tuple(int(ZZ(x) % MODULUS) for x in v)


def in_e8(v):
    if len(v) != E8_BLOCK_DIM:
        return False
    if not is_integral_vector(v):
        return False
    return parity_tuple(v) in RM13_CODEWORD_SET


def in_standard_half_integer_e8(v):
    # Reference model for the usual E8 description:
    #   {x in Z^8 : sum x_i even} union
    #   {x in (Z + 1/2)^8 : sum x_i even}.
    # The sampler uses the Construction-A model above, but this helper makes
    # the familiar half integer structure explicit for checks and exposition.
    if len(v) != E8_BLOCK_DIM:
        return False
    q = [QQ(x) for x in v]
    if all(x in ZZ for x in q):
        return ZZ(sum(q)) % 2 == 0
    doubled = [2 * x for x in q]
    if all(x in ZZ and ZZ(x) % 2 != 0 for x in doubled):
        total = sum(q)
        return total in ZZ and ZZ(total) % 2 == 0
    return False


def split_e8_blocks(v, k=None):
    n = len(v)
    if n % E8_BLOCK_DIM != 0:
        raise ValueError("vector length must be a multiple of 8")
    inferred_k = ZZ(n // E8_BLOCK_DIM)
    if k is not None and ZZ(k) != inferred_k:
        raise ValueError("vector length does not match k")
    return [
        vector(QQ, list(v[E8_BLOCK_DIM * i:E8_BLOCK_DIM * (i + 1)]))
        for i in range(inferred_k)
    ]


def in_e8_power(v, k=None):
    return all(in_e8(block) for block in split_e8_blocks(v, k))


def in_2e8(v):
    if len(v) != E8_BLOCK_DIM:
        return False
    if not is_integral_vector(v):
        return False
    if any(ZZ(x) % 2 != 0 for x in v):
        return False
    return in_e8(vector(ZZ, [ZZ(x) // 2 for x in v]))


def in_2e8_power(v, k=None):
    return all(in_2e8(block) for block in split_e8_blocks(v, k))


def in_coset_mod_2e8(v, rep, k=None):
    if len(v) != len(rep):
        return False
    diff = vector(QQ, [QQ(v[i]) - QQ(rep[i]) for i in range(len(v))])
    return in_2e8_power(diff, k)


def squared_norm(v):
    return sum(QQ(x) * QQ(x) for x in v)


def half_shift(rep):
    return vector(QQ, [QQ(x) / 2 for x in rep])


def e8_mod_2e8_representatives():
    reps = []
    for bits in itertools.product([0, 1], repeat=E8_BLOCK_DIM):
        reps.append(coset_rep_from_bits(bits, k=1))
    return reps


def coset_rep_from_bits(bits, k=DEFAULT_BLOCKS):
    bits = [ZZ(b) for b in bits]
    k = ZZ(k)
    if len(bits) != challenge_bit_count(k):
        raise ValueError("expected exactly 8k bits")
    if any(b not in (0, 1) for b in bits):
        raise ValueError("coset representative bits must be binary")

    out = [ZZ(0)] * e8_dimension(k)
    for block in range(k):
        r = vector(ZZ, [0] * E8_BLOCK_DIM)
        for i in range(E8_BLOCK_DIM):
            if bits[E8_BLOCK_DIM * block + i]:
                r += E8_BASIS.row(i)
        for i in range(E8_BLOCK_DIM):
            out[E8_BLOCK_DIM * block + i] = ZZ(r[i])
    return vector(ZZ, out)


def direct_sum_e8_basis(k=DEFAULT_BLOCKS):
    k = ZZ(k)
    if k <= 0:
        raise ValueError("k must be positive")
    dim = e8_dimension(k)
    basis = matrix(ZZ, dim, dim)
    for block in range(k):
        offset = E8_BLOCK_DIM * block
        for i in range(E8_BLOCK_DIM):
            for j in range(E8_BLOCK_DIM):
                basis[offset + i, offset + j] = E8_BASIS[i, j]
    return basis


def serialize_integer_matrix(mat):
    rows = []
    for i in range(mat.nrows()):
        rows.append(",".join(str(ZZ(mat[i, j])) for j in range(mat.ncols())))
    return (";".join(rows)).encode("ascii")
