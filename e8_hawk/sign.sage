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
load(os.path.join(E8_HAWK_DIR, "hash_to_target.sage"))
load(os.path.join(E8_HAWK_DIR, "symbreak.sage"))
load(os.path.join(E8_HAWK_DIR, "material.sage"))

from sage.all import QQ, ZZ, vector


def signature_from_h_and_w(h, w, k=DEFAULT_BLOCKS):
    if len(h) != len(w):
        raise ValueError("h and w must have the same length")
    if not in_coset_mod_2e8(w, h, k=k):
        raise ValueError("w must be congruent to h modulo 2E8^k")
    s = vector(QQ, [(QQ(h[i]) - QQ(w[i])) / 2 for i in range(len(h))])
    if not in_e8_power(s, k):
        raise ValueError("derived signature vector is not in E8^k")
    return s


def compact_signature(signature):
    return {
        "salt": signature["salt"],
        "s": signature["s"],
    }


def sign(
    message,
    secret_material,
    salt=None,
    norm_bound=None,
    max_attempts=MAX_REJECTION_ATTEMPTS,
):
    validate_secret_material(secret_material)
    public_material = secret_material_public(secret_material)
    validate_public_material(public_material)
    k = public_material_k(public_material)
    sigma_sign = secret_material_sigma_sign(secret_material)
    sigma_verify = public_material_sigma_verify(public_material)
    public_context = public_material_context(public_material)

    target = hash_to_target(message, salt=salt, k=k, public_context=public_context)
    h = target["h"]
    bound = public_material_norm_bound(public_material, k, sigma_verify) if norm_bound is None else QQ(norm_bound)

    attempts = 0
    norm_rejections = 0
    zero_rejections = 0
    rejected_norms = []

    while True:
        attempts += 1
        if max_attempts is not None and attempts > max_attempts:
            raise RuntimeError("rejection loop exceeded max_attempts")

        sampled_w = sample_w_from_secret_material(secret_material, h, k, sigma_sign)
        if not public_material_witness_membership(public_material, sampled_w):
            raise RuntimeError("sampler returned a witness outside the public lattice")

        sampled_norm_sq = public_material_norm_sq(public_material, sampled_w)
        if sampled_norm_sq > bound:
            norm_rejections += 1
            rejected_norms.append(sampled_norm_sq)
            continue

        try:
            w, symbreak_flipped = canonical_representative(sampled_w)
        except ValueError:
            zero_rejections += 1
            continue

        if not sym_break(w):
            raise RuntimeError("canonical representative failed sym-break")
        if not public_material_witness_membership(public_material, w):
            raise RuntimeError("canonical representative is outside the public lattice")
        if public_material_norm_sq(public_material, w) > bound:
            raise RuntimeError("canonical representative exceeded the norm bound")

        s = signature_from_h_and_w(h, w, k=k)
        if not public_material_signature_membership(public_material, s):
            raise RuntimeError("derived signature vector is outside the signature module")
        return {
            "message": message,
            "salt": target["salt"],
            "public_context": target["public_context"],
            "k": ZZ(k),
            "sigma_sign": QQ(sigma_sign),
            "sigma_verify": QQ(sigma_verify),
            "norm_bound_sq": bound,
            "public_material": public_material,
            "h": h,
            "h_bits": target["bits"],
            "sampled_w": sampled_w,
            "w": w,
            "s": s,
            "attempts": attempts,
            "norm_rejections": norm_rejections,
            "zero_rejections": zero_rejections,
            "rejected_norms": rejected_norms,
            "symbreak_flipped": symbreak_flipped,
            "w_in_public_lattice": public_material_witness_membership(public_material, w),
            "w_norm_sq": public_material_norm_sq(public_material, w),
            "norm_within_bound": public_material_norm_sq(public_material, w) <= bound,
            "symbreak": sym_break(w),
            "relation_h_minus_2s": h - 2 * s == w,
        }
