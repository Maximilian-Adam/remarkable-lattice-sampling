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
load(os.path.join(E8_HAWK_DIR, "e8_lattice.sage"))

from sage.all import QQ, ZZ, vector
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler as DGaussZ
from math import isfinite
from mpmath import mp, exp, jtheta, log, pi, sqrt
from random import SystemRandom

mp.dps = MP_DPS

_DEFAULT_SAMPLER_RNG = SystemRandom()


def default_sampler_rng():
    return _DEFAULT_SAMPLER_RNG


def _rng_random(rng):
    rng = default_sampler_rng() if rng is None else rng
    if isinstance(rng, dict) and "random" in rng:
        value = rng["random"]()
    elif hasattr(rng, "random"):
        value = rng.random()
    elif callable(rng):
        value = rng()
    else:
        raise ValueError("rng must provide random() or be a zero-argument callback")
    value = float(value)
    if value < 0.0 or value >= 1.0:
        raise ValueError("rng random value must lie in [0, 1)")
    return value


def _rng_coord_sampler(rng):
    if isinstance(rng, dict):
        return rng.get("sample_coord_2z_coset", rng.get("sample_coord", None))
    if hasattr(rng, "sample_coord_2z_coset"):
        return rng.sample_coord_2z_coset
    if hasattr(rng, "sample_coord"):
        return rng.sample_coord
    return None


def nome(s):
    return exp(-pi / (float(s) ** 2))


def rho_2z_shift(alpha, s):
    q = nome(s)
    alpha_f = float(alpha)
    z = -2j * pi * alpha_f / (float(s) ** 2)
    value = (q ** (alpha_f ** 2)) * jtheta(3, z, q ** 4)
    return float(value.real)


def codeword_weights_shifted_e8(s, shift):
    shift_f = [float(a) for a in shift]
    log_weights = []
    for c in RM13_CODEWORDS:
        lw = 0.0
        for i in range(E8_BLOCK_DIM):
            lw += float(log(rho_2z_shift(c[i] + shift_f[i], s)))
        log_weights.append(lw)
    m = max(log_weights)
    return [float(exp(w - m)) for w in log_weights]


def sample_index_from_weights(weights, rng=None):
    weights = [float(w) for w in weights]
    if not weights:
        raise ValueError("weights must be non-empty")
    if any((not isfinite(w)) or w < 0.0 for w in weights):
        raise ValueError("weights must be finite and non-negative")
    total = sum(weights)
    if total <= 0.0:
        raise ValueError("at least one weight must be positive")
    threshold = _rng_random(rng) * total
    cumulative = 0.0
    for i, w in enumerate(weights):
        cumulative += w
        if threshold < cumulative:
            return i
    return len(weights) - 1


def sample_coord_2z_coset(alpha, s, rng=None):
    coord_sampler = _rng_coord_sampler(rng)
    if coord_sampler is not None:
        return QQ(coord_sampler(alpha, s))
    # Sage's integer sampler uses exp(-(x-c)^2/(2*sigma^2)); this converts
    # from rho_s(x)=exp(-pi*x^2/s^2) on the 2Z+alpha coset.
    sigma = float(s) / (2 * float(sqrt(2 * pi)))
    center = -QQ(alpha) / 2
    z = DGaussZ(sigma=sigma, c=center)()
    return ZZ(2 * z) + QQ(alpha)


def sample_shifted_e8(s, shift, rng=None):
    if len(shift) != E8_BLOCK_DIM:
        raise ValueError("shift must have length 8")
    weights = codeword_weights_shifted_e8(s, shift)
    idx = sample_index_from_weights(weights, rng=rng)
    c = RM13_CODEWORDS[idx]
    x = [sample_coord_2z_coset(c[i] + shift[i], s, rng=rng) for i in range(E8_BLOCK_DIM)]
    return vector(QQ, x), c


def sample_coset_rep_plus_2e8(width, rep, rng=None):
    if not in_e8(rep):
        raise ValueError("rep must lie in E8")
    y, _ = sample_shifted_e8(QQ(width) / 2, half_shift(rep), rng=rng)
    out = vector(QQ, [2 * yi for yi in y])
    if not in_coset_mod_2e8(out, rep, k=1):
        raise RuntimeError("internal sampler produced a vector outside rep + 2E8")
    return out


def sample_e8_power_coset(width, rep, k=None, rng=None):
    blocks = split_e8_blocks(rep, k)
    samples = [sample_coset_rep_plus_2e8(width, block, rng=rng) for block in blocks]
    out = []
    for sample in samples:
        out.extend(sample)
    return vector(QQ, out)


def sample_candidate_w(h, k=DEFAULT_BLOCKS, sigma_sign=SIGMA_SIGN, rng=None):
    # HAWK-style stage: draw the public witness from h + 2E8^k at width
    # 2*sigma_sign, then let sign.sage apply norm rejection and sym-break.
    width = 2 * QQ(sigma_sign)
    w = sample_e8_power_coset(width, h, k=k, rng=rng)
    if not in_e8_power(w, k):
        raise RuntimeError("candidate is not in E8^k")
    if not in_coset_mod_2e8(w, h, k=k):
        raise RuntimeError("candidate is not congruent to h modulo 2E8^k")
    return w
