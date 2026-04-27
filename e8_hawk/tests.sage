import os
import sys

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
load(os.path.join(E8_HAWK_DIR, "sampler.sage"))
load(os.path.join(E8_HAWK_DIR, "symbreak.sage"))
load(os.path.join(E8_HAWK_DIR, "material.sage"))
load(os.path.join(E8_HAWK_DIR, "sign.sage"))
load(os.path.join(E8_HAWK_DIR, "verify.sage"))

from sage.all import QQ, ZZ, identity_matrix, vector


def assert_true(condition, label):
    if not condition:
        raise AssertionError(label)


def make_signature_from_w(message, salt, w, public_material, norm_bound=None):
    k = public_material["k"]
    target = hash_to_target(message, salt=salt, k=k, public_context=public_material["public_context"])
    h = target["h"]
    s = signature_from_h_and_w(h, w, k=k)
    return {
        "salt": target["salt"],
        "s": s,
        "norm_bound_sq": public_material["norm_bound_sq"] if norm_bound is None else QQ(norm_bound),
    }


def test_membership_checks():
    zero = vector(ZZ, [0] * 8)
    assert_true(in_e8(zero), "zero should be in E8")
    assert_true(in_e8(E8_BASIS.row(0)), "basis row should be in E8")
    assert_true(not in_e8(vector(ZZ, [1, 0, 0, 0, 0, 0, 0, 0])), "unit vector should not be in this E8 model")
    assert_true(in_e8_power(vector(ZZ, list(zero) + list(E8_BASIS.row(1))), k=2), "blockwise E8^2 membership failed")
    assert_true(in_standard_half_integer_e8(vector(QQ, [QQ(1) / 2] * 8)), "standard half-integer E8 check failed")


def test_hashing_consistency():
    salt = bytes.fromhex("00112233445566778899aabbccddeeff")
    public_material = make_toy_e8_public_material(k=2)
    a = hash_to_target("message", salt=salt, k=2, public_context=public_material["public_context"])
    b = hash_to_target("message", salt=salt, k=2, public_context=public_material["public_context"])
    c = hash_to_target("message!", salt=salt, k=2, public_context=public_material["public_context"])
    d = hash_to_target("message", salt=salt, k=2, public_context=b"different-public-context")
    assert_true(a["h"] == b["h"], "hash_to_target must be deterministic")
    assert_true(a["bits"] == b["bits"], "hash bits must be deterministic")
    assert_true(a["h"] != c["h"], "different messages should change the target in this test")
    assert_true(a["h"] != d["h"], "different public context should change the target in this test")
    assert_true(in_e8_power(a["h"], k=2), "hash target should lie in E8^k")


def test_valid_signature_verifies():
    public_material = make_toy_e8_public_material(k=1)
    secret_material = make_toy_e8_secret_material(public_material)
    sig = sign("valid signature", secret_material, salt=bytes.fromhex("10112233445566778899aabbccddeeff"))
    compact = compact_signature(sig)
    assert_true(verify("valid signature", compact, public_material), "valid signature did not verify")
    details = verify("valid signature", compact, public_material, return_details=True)
    assert_true(details["checks"]["w_in_public_lattice"], "verified witness should be in E8^k")
    assert_true(details["checks"]["norm_bound"], "verified witness should satisfy norm bound")
    assert_true(details["checks"]["symbreak"], "verified witness should satisfy sym-break")
    assert_true(sig["relation_h_minus_2s"], "signing relation w = h - 2s failed")


def test_q_public_norm_path():
    message = "q-norm"
    salt = bytes.fromhex("12112233445566778899aabbccddeeff")
    k = ZZ(1)
    q_scale = QQ(3)
    q = q_scale * identity_matrix(QQ, e8_dimension(k))
    public_material = make_toy_e8_public_material(k=k, Q=q)
    assert_true("norm_sq" not in public_material, "toy material should use Q path for this test")

    target = hash_to_target(message, salt=salt, k=k, public_context=public_material["public_context"])
    w = target["h"] + 2 * E8_BASIS.row(0)
    w, _ = canonical_representative(w)
    exact_q_norm = public_material_norm_sq(public_material, w)
    assert_true(exact_q_norm == q_scale * squared_norm(w), "Q norm path did not scale the norm")

    exact_public_material = make_toy_e8_public_material(k=k, Q=q, norm_bound=exact_q_norm)
    sig = make_signature_from_w(message, salt, w, exact_public_material, norm_bound=exact_q_norm)
    assert_true(verify(message, sig, exact_public_material), "witness at exact Q norm bound should verify")

    tight_public_material = make_toy_e8_public_material(k=k, Q=q, norm_bound=exact_q_norm - 1)
    assert_true(not verify(message, sig, tight_public_material), "witness above a tightened Q bound should fail")


def test_mismatched_public_material_fails():
    public_material = make_toy_e8_public_material(k=1)
    secret_material = make_toy_e8_secret_material(public_material)
    sig = compact_signature(sign("material mismatch", secret_material, salt=bytes.fromhex("13112233445566778899aabbccddeeff")))

    mismatched_q = 1000 * identity_matrix(QQ, e8_dimension(1))
    mismatched_public_material = make_toy_e8_public_material(
        k=1,
        Q=mismatched_q,
        public_context=b"different-public-material-v1",
    )
    details = verify("material mismatch", sig, mismatched_public_material, return_details=True)
    assert_true(not details["valid"], "signature should fail under mismatched public material")
    assert_true(not details["checks"]["norm_bound"], "mismatched Q should fail the public norm bound")


def test_material_is_required():
    public_material = make_toy_e8_public_material(k=1)
    secret_material = make_toy_e8_secret_material(public_material)
    raised = False
    try:
        sign("missing secret material")
    except TypeError:
        raised = True
    except ValueError:
        raised = True
    assert_true(raised, "signing should require explicit secret material")

    sig = compact_signature(sign("missing public material", secret_material, salt=bytes.fromhex("11112233445566778899aabbccddeeff")))
    raised = False
    try:
        verify("missing public material", sig)
    except TypeError:
        raised = True
    except ValueError:
        raised = True
    assert_true(raised, "verification should require explicit public material")


def test_tampered_signature_fails():
    public_material = make_toy_e8_public_material(k=1)
    secret_material = make_toy_e8_secret_material(public_material)
    sig = compact_signature(sign("tamper", secret_material, salt=bytes.fromhex("20112233445566778899aabbccddeeff")))
    tampered = dict(sig)
    tampered["s"] = vector(QQ, sig["s"]) + vector(QQ, [QQ(1) / 2, 0, 0, 0, 0, 0, 0, 0])
    details = verify("tamper", tampered, public_material, return_details=True)
    assert_true(not details["valid"], "tampered signature vector should fail verification")


def test_negated_signature_fails_by_symbreak():
    public_material = make_toy_e8_public_material(k=1)
    secret_material = make_toy_e8_secret_material(public_material)
    sig_full = sign("symbreak", secret_material, salt=bytes.fromhex("30112233445566778899aabbccddeeff"))
    h = sig_full["h"]
    s = sig_full["s"]
    transformed = compact_signature(sig_full)
    transformed["s"] = h - s
    details = verify("symbreak", transformed, public_material, return_details=True)
    assert_true(not details["valid"], "transformed signature should fail")
    assert_true(not details["checks"]["symbreak"], "transformed signature should fail specifically by sym-break")


def test_oversized_witness_fails():
    message = "oversized"
    salt = bytes.fromhex("40112233445566778899aabbccddeeff")
    public_material = make_toy_e8_public_material(k=1)
    target = hash_to_target(message, salt=salt, k=1, public_context=public_material["public_context"])
    huge_direction = 100 * E8_BASIS.row(0)
    w = target["h"] + 2 * huge_direction
    w, _ = canonical_representative(w)
    sig = make_signature_from_w(message, salt, w, public_material)
    details = verify(message, sig, public_material, return_details=True)
    assert_true(not details["valid"], "oversized witness should fail")
    assert_true(details["checks"]["w_in_public_lattice"], "oversized witness should still be in E8")
    assert_true(not details["checks"]["norm_bound"], "oversized witness should fail the norm bound")


def test_non_e8_witness_fails():
    message = "non-e8"
    salt = bytes.fromhex("50112233445566778899aabbccddeeff")
    public_material = make_toy_e8_public_material(k=1)
    target = hash_to_target(message, salt=salt, k=1, public_context=public_material["public_context"])
    bad_w = vector(QQ, target["h"]) + vector(QQ, [1, 0, 0, 0, 0, 0, 0, 0])
    sig = {
        "salt": target["salt"],
        "s": vector(QQ, [(QQ(target["h"][i]) - bad_w[i]) / 2 for i in range(8)]),
    }
    details = verify(message, sig, public_material, return_details=True)
    assert_true(not details["valid"], "non-E8 witness should fail")
    assert_true(not details["checks"]["w_in_public_lattice"], "non-E8 witness check did not fail")


def test_norm_boundary_cases():
    message = "boundary"
    salt = bytes.fromhex("60112233445566778899aabbccddeeff")
    public_material = make_toy_e8_public_material(k=1)
    target = hash_to_target(message, salt=salt, k=1, public_context=public_material["public_context"])
    w = target["h"] + 2 * E8_BASIS.row(0)
    w, _ = canonical_representative(w)
    exact = squared_norm(w)
    exact_public_material = make_toy_e8_public_material(k=1, norm_bound=exact)
    sig = make_signature_from_w(message, salt, w, exact_public_material, norm_bound=exact)
    assert_true(verify(message, sig, exact_public_material), "witness at exact norm bound should verify")
    tight_public_material = make_toy_e8_public_material(k=1, norm_bound=exact - 1)
    assert_true(not verify(message, sig, tight_public_material), "witness just above a tightened bound should fail")


def run_all_tests():
    tests = [
        test_membership_checks,
        test_hashing_consistency,
        test_valid_signature_verifies,
        test_q_public_norm_path,
        test_mismatched_public_material_fails,
        test_material_is_required,
        test_tampered_signature_fails,
        test_negated_signature_fails_by_symbreak,
        test_oversized_witness_fails,
        test_non_e8_witness_fails,
        test_norm_boundary_cases,
    ]
    for test in tests:
        test()
        print("PASS", test.__name__)
    print("all E8-HAWK signing-stage tests passed")


run_all_tests()
