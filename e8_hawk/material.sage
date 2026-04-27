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
load(os.path.join(E8_HAWK_DIR, "sampler.sage"))

from sage.all import QQ, ZZ, identity_matrix, matrix, vector


def public_material_k(public_material):
    return ZZ(public_material.get("k", DEFAULT_BLOCKS))


def public_material_context(public_material):
    return public_material.get("public_context", b"")


def public_material_sigma_verify(public_material):
    return QQ(public_material.get("sigma_verify", SIGMA_VERIFY))


def public_material_norm_bound(public_material, k=None, sigma_verify=None):
    if "norm_bound_sq" in public_material:
        return QQ(public_material["norm_bound_sq"])
    if k is None:
        k = public_material_k(public_material)
    if sigma_verify is None:
        sigma_verify = public_material_sigma_verify(public_material)
    return norm_bound_sq(k, sigma_verify)


def public_material_norm_sq(public_material, w):
    if "norm_sq" in public_material:
        return QQ(public_material["norm_sq"](w, public_material))
    if "Q" in public_material:
        w_q = vector(QQ, w)
        return QQ(w_q.dot_product(w_q * public_material["Q"]))
    return squared_norm(w)


def public_material_witness_membership(public_material, w):
    if "witness_membership" in public_material:
        return bool(public_material["witness_membership"](w, public_material))
    return in_e8_power(w, public_material_k(public_material))


def public_material_signature_membership(public_material, s):
    if "signature_membership" in public_material:
        return bool(public_material["signature_membership"](s, public_material))
    return in_e8_power(s, public_material_k(public_material))


def validate_public_material(public_material):
    if public_material is None:
        raise ValueError("public_material is required")
    k = public_material_k(public_material)
    if k <= 0:
        raise ValueError("public_material k must be positive")
    if "Q" in public_material:
        q = public_material["Q"]
        dim = e8_dimension(k)
        if q.nrows() != dim or q.ncols() != dim:
            raise ValueError("public_material Q must be an 8k x 8k matrix")
    public_material_norm_bound(public_material, k, public_material_sigma_verify(public_material))
    public_material_context(public_material)
    return True


def secret_material_public(secret_material):
    if secret_material is None:
        raise ValueError("sign() requires secret_material; signing-stage code must not generate keys")
    if "public_material" in secret_material:
        return secret_material["public_material"]
    if "public" in secret_material:
        return secret_material["public"]
    raise ValueError("secret_material must contain public_material")


def secret_material_sigma_sign(secret_material):
    return QQ(secret_material.get("sigma_sign", SIGMA_SIGN))


def sample_w_from_secret_material(secret_material, h, k, sigma_sign):
    sampler = secret_material.get("sample_w", secret_material.get("sampler", None))
    if sampler is None:
        raise ValueError("secret_material must provide sample_w(h, k, sigma_sign)")
    return sampler(h, k=k, sigma_sign=sigma_sign)


def validate_secret_material(secret_material):
    public_material = secret_material_public(secret_material)
    validate_public_material(public_material)
    if "sample_w" not in secret_material and "sampler" not in secret_material:
        raise ValueError("secret_material must provide sample_w or sampler")
    secret_material_sigma_sign(secret_material)
    return True


def toy_e8_public_context(k, label=b"toy-e8-public-material-v1"):
    return label + int(ZZ(k)).to_bytes(4, "big")


def toy_e8_witness_membership(w, public_material):
    return in_e8_power(w, public_material_k(public_material))


def toy_e8_signature_membership(s, public_material):
    return in_e8_power(s, public_material_k(public_material))


def toy_e8_sample_w(h, k, sigma_sign):
    return sample_candidate_w(h, k=k, sigma_sign=sigma_sign)


def make_toy_e8_public_material(
    k=DEFAULT_BLOCKS,
    sigma_verify=SIGMA_VERIFY,
    norm_bound=None,
    public_context=None,
    Q=None,
):
    # This is an adapter fixture, not HAWK KeyGen. It exposes the
    # public data the signing stage needs: dimensions, transcript context,
    # public norm object Q, and membership predicates.
    k = ZZ(k)
    if public_context is None:
        public_context = toy_e8_public_context(k)
    if Q is None:
        Q = identity_matrix(QQ, e8_dimension(k))
    public_material = {
        "name": "toy-e8-public-material",
        "k": k,
        "sigma_verify": QQ(sigma_verify),
        "norm_bound_sq": norm_bound_sq(k, sigma_verify) if norm_bound is None else QQ(norm_bound),
        "public_context": public_context,
        "Q": matrix(QQ, Q),
        "witness_membership": toy_e8_witness_membership,
        "signature_membership": toy_e8_signature_membership,
    }
    validate_public_material(public_material)
    return public_material


def make_toy_e8_secret_material(public_material, sigma_sign=SIGMA_SIGN):
    # This is the toy E8 sampling adapter used by tests and experiments. A
    # real HAWK integration needs to provide its own external secret material.
    secret_material = {
        "name": "toy-e8-secret-material",
        "public_material": public_material,
        "sigma_sign": QQ(sigma_sign),
        "sample_w": toy_e8_sample_w,
    }
    validate_secret_material(secret_material)
    return secret_material
