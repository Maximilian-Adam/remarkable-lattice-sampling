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


def verify(
    message,
    signature,
    public_material,
    norm_bound=None,
    return_details=False,
):
    validate_public_material(public_material)

    k = public_material_k(public_material)
    sigma_verify = public_material_sigma_verify(public_material)
    public_context = public_material_context(public_material)
    bound = public_material_norm_bound(public_material, k, sigma_verify) if norm_bound is None else QQ(norm_bound)
    target = hash_to_target(
        message,
        salt=signature["salt"],
        k=k,
        public_context=public_context,
    )
    h = target["h"]

    try:
        s = vector(QQ, list(signature["s"]))
    except Exception:
        details = {"valid": False, "reason": "signature vector is not a vector"}
        return details if return_details else False

    checks = {}
    checks["length"] = len(s) == e8_dimension(k)
    if not checks["length"]:
        details = {"valid": False, "reason": "signature vector has wrong length", "checks": checks}
        return details if return_details else False

    w = h - 2 * s
    checks["s_in_public_signature_module"] = public_material_signature_membership(public_material, s)
    checks["w_in_public_lattice"] = public_material_witness_membership(public_material, w)
    checks["same_2e8_coset"] = in_coset_mod_2e8(w, h, k=k)
    checks["symbreak"] = sym_break(w) if checks["w_in_public_lattice"] else False
    w_norm_sq = public_material_norm_sq(public_material, w)
    checks["norm_bound"] = w_norm_sq <= bound

    valid = all(checks.values())
    details = {
        "valid": valid,
        "checks": checks,
        "h": h,
        "w": w,
        "s": s,
        "w_norm_sq": w_norm_sq,
        "norm_bound_sq": bound,
        "salt": target["salt"],
        "k": ZZ(k),
        "sigma_verify": QQ(sigma_verify),
    }
    return details if return_details else valid
