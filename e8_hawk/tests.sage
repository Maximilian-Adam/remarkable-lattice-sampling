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


def fixed_salt(tag, salt_bytes=SALT_BYTES):
    return bytes((int(tag) + i) % 256 for i in range(int(salt_bytes)))


def make_signature_from_w(message, salt, w, public_material, norm_bound=None):
    k = public_material["k"]
    target = hash_to_target(
        message,
        salt=salt,
        k=k,
        public_context=public_material["public_context"],
        salt_bytes=public_material_salt_bytes(public_material),
    )
    h = target["h"]
    s = signature_from_h_and_w(h, w, k=k, public_material=public_material)
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
    salt = fixed_salt(0)
    public_material = make_toy_e8_public_material(k=2)
    salt_bytes = public_material_salt_bytes(public_material)
    a = hash_to_target("message", salt=salt, k=2, public_context=public_material["public_context"], salt_bytes=salt_bytes)
    b = hash_to_target("message", salt=salt, k=2, public_context=public_material["public_context"], salt_bytes=salt_bytes)
    c = hash_to_target("message!", salt=salt, k=2, public_context=public_material["public_context"], salt_bytes=salt_bytes)
    d = hash_to_target("message", salt=salt, k=2, public_context=b"different-public-context", salt_bytes=salt_bytes)
    assert_true(a["h"] == b["h"], "hash_to_target must be deterministic")
    assert_true(a["bits"] == b["bits"], "hash bits must be deterministic")
    assert_true(a["h"] != c["h"], "different messages should change the target in this test")
    assert_true(a["h"] != d["h"], "different public context should change the target in this test")
    assert_true(in_e8_power(a["h"], k=2), "hash target should lie in E8^k")


def test_pinned_basis_and_representative_order():
    expected_basis = (
        "1,0,0,1,0,1,1,0;"
        "0,1,0,1,0,1,0,1;"
        "0,0,1,1,0,0,1,1;"
        "0,0,0,2,0,0,0,0;"
        "0,0,0,0,1,1,1,1;"
        "0,0,0,0,0,2,0,0;"
        "0,0,0,0,0,0,2,0;"
        "0,0,0,0,0,0,0,2"
    )
    assert_true(serialize_integer_matrix(E8_BASIS).decode("ascii") == expected_basis, "E8 basis order changed")
    bits = [1, 0, 1, 0, 1, 0, 1, 0]
    rep = coset_rep_from_bits(bits, k=1)
    assert_true(list(rep) == [1, 0, 1, 2, 1, 2, 5, 2], "bit-to-representative order changed")


def test_salt_length_is_fixed():
    valid = fixed_salt(1)
    assert_true(len(normalize_salt(None)) == SALT_BYTES, "generated salt has wrong length")
    assert_true(normalize_salt(valid.hex()) == valid, "valid hex salt failed")
    for bad in (b"short", "00", bytes(range(SALT_BYTES + 1))):
        raised = False
        try:
            normalize_salt(bad)
        except ValueError:
            raised = True
        assert_true(raised, "salt with wrong length should fail")


def test_hawk_parameter_sets():
    p256 = hawk_parameter_set("HAWK-256")
    p512 = hawk_parameter_set("HAWK-512")
    p1024 = hawk_parameter_set("HAWK-1024")
    assert_true(p256["degree"] == 256 and p256["e8_blocks"] == 64, "HAWK-256 dimension mapping is wrong")
    assert_true(p512["degree"] == 512 and p512["e8_blocks"] == 128, "HAWK-512 dimension mapping is wrong")
    assert_true(p1024["degree"] == 1024 and p1024["e8_blocks"] == 256, "HAWK-1024 dimension mapping is wrong")
    assert_true(p256["salt_bytes"] == 14 and p512["salt_bytes"] == 24 and p1024["salt_bytes"] == 40, "HAWK salt lengths are wrong")
    assert_true(p512["sigma_sign"] == QQ(1278) / QQ(1000), "HAWK-512 sigma_sign is wrong")
    assert_true(p512["sigma_verify"] == QQ(1425) / QQ(1000), "HAWK-512 sigma_verify is wrong")
    assert_true(DEFAULT_BLOCKS == HAWK512_BLOCKS and SALT_BYTES == 24, "default profile should be HAWK-512")

    summary = parameter_summary(parameter_set="HAWK-1024")
    assert_true(summary["k"] == 256 and summary["dimension"] == 2048, "HAWK-1024 summary dimension is wrong")
    assert_true(summary["salt_bytes"] == 40, "HAWK-1024 summary salt length is wrong")

    public_material = make_toy_e8_public_material(parameter_set="HAWK-256")
    secret_material = make_toy_e8_secret_material(public_material)
    assert_true(public_material["k"] == 64, "HAWK-256 toy material should use 64 E8 blocks")
    assert_true(public_material_salt_bytes(public_material) == 14, "HAWK-256 toy material salt length is wrong")
    assert_true(public_material_sigma_verify(public_material) == QQ(1042) / QQ(1000), "HAWK-256 verifier width is wrong")
    assert_true(secret_material_sigma_sign(secret_material) == QQ(1010) / QQ(1000), "HAWK-256 signer width is wrong")


class FixedRandom:
    def __init__(self, values):
        self.values = list(values)

    def random(self):
        if not self.values:
            raise RuntimeError("fixed random stream exhausted")
        return self.values.pop(0)


def test_sampler_rng_hooks():
    assert_true(sample_index_from_weights([1.0, 3.0], rng=FixedRandom([0.0])) == 0, "weighted index lower bucket failed")
    assert_true(sample_index_from_weights([1.0, 3.0], rng=FixedRandom([0.99])) == 1, "weighted index upper bucket failed")

    calls = []
    def coord(alpha, width):
        calls.append((QQ(alpha), QQ(width)))
        return QQ(alpha)

    rng = {"random": FixedRandom([0.0]).random, "sample_coord_2z_coset": coord}
    sample, codeword = sample_shifted_e8(QQ(1), vector(QQ, [0] * E8_BLOCK_DIM), rng=rng)
    assert_true(codeword == RM13_CODEWORDS[0], "rng hook should select the first codeword at random value 0")
    assert_true(sample == vector(QQ, [0] * E8_BLOCK_DIM), "coordinate sampler hook was not used")
    assert_true(len(calls) == E8_BLOCK_DIM, "coordinate sampler hook should be called once per coordinate")


def test_valid_signature_verifies():
    public_material = make_toy_e8_public_material(k=1)
    secret_material = make_toy_e8_secret_material(public_material)
    sig = sign("valid signature", secret_material, salt=fixed_salt(0x10))
    compact = compact_signature(sig)
    assert_true(verify("valid signature", compact, public_material), "valid signature did not verify")
    details = verify("valid signature", compact, public_material, return_details=True)
    assert_true(details["checks"]["w_in_public_lattice"], "verified witness should be in E8^k")
    assert_true(details["checks"]["norm_bound"], "verified witness should satisfy norm bound")
    assert_true(details["checks"]["symbreak"], "verified witness should satisfy sym-break")
    assert_true(sig["relation_h_minus_2s"], "signing relation w = h - 2s failed")


def test_q_public_norm_path():
    message = "q-norm"
    salt = fixed_salt(0x12)
    k = ZZ(1)
    q_scale = QQ(3)
    q = q_scale * identity_matrix(QQ, e8_dimension(k))
    public_material = make_toy_e8_public_material(k=k, Q=q)
    assert_true("norm_sq" not in public_material, "toy material should use Q path for this test")

    target = hash_to_target(message, salt=salt, k=k, public_context=public_material["public_context"], salt_bytes=public_material_salt_bytes(public_material))
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
    sig = compact_signature(sign("material mismatch", secret_material, salt=fixed_salt(0x13)))

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

    sig = compact_signature(sign("missing public material", secret_material, salt=fixed_salt(0x11)))
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
    sig = compact_signature(sign("tamper", secret_material, salt=fixed_salt(0x20)))
    tampered = dict(sig)
    tampered["s"] = vector(QQ, sig["s"]) + vector(QQ, [QQ(1) / 2, 0, 0, 0, 0, 0, 0, 0])
    details = verify("tamper", tampered, public_material, return_details=True)
    assert_true(not details["valid"], "tampered signature vector should fail verification")


def test_negated_signature_fails_by_symbreak():
    public_material = make_toy_e8_public_material(k=1)
    secret_material = make_toy_e8_secret_material(public_material)
    sig_full = sign("symbreak", secret_material, salt=fixed_salt(0x30))
    h = sig_full["h"]
    s = sig_full["s"]
    transformed = compact_signature(sig_full)
    transformed["s"] = h - s
    details = verify("symbreak", transformed, public_material, return_details=True)
    assert_true(not details["valid"], "transformed signature should fail")
    assert_true(not details["checks"]["symbreak"], "transformed signature should fail specifically by sym-break")


def test_oversized_witness_fails():
    message = "oversized"
    salt = fixed_salt(0x40)
    public_material = make_toy_e8_public_material(k=1)
    target = hash_to_target(message, salt=salt, k=1, public_context=public_material["public_context"], salt_bytes=public_material_salt_bytes(public_material))
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
    salt = fixed_salt(0x50)
    public_material = make_toy_e8_public_material(k=1)
    target = hash_to_target(message, salt=salt, k=1, public_context=public_material["public_context"], salt_bytes=public_material_salt_bytes(public_material))
    bad_w = vector(QQ, target["h"]) + vector(QQ, [1, 0, 0, 0, 0, 0, 0, 0])
    sig = {
        "salt": target["salt"],
        "s": vector(QQ, [(QQ(target["h"][i]) - bad_w[i]) / 2 for i in range(8)]),
    }
    details = verify(message, sig, public_material, return_details=True)
    assert_true(not details["valid"], "non-E8 witness should fail")
    assert_true(not details["checks"]["w_in_public_lattice"], "non-E8 witness check did not fail")


def test_public_coset_relation_is_material_controlled():
    message = "public-coset"
    salt = fixed_salt(0x51)
    public_material = make_toy_e8_public_material(k=1)
    target = hash_to_target(message, salt=salt, k=1, public_context=public_material["public_context"], salt_bytes=public_material_salt_bytes(public_material))
    w = target["h"] + 2 * E8_BASIS.row(0)
    w, _ = canonical_representative(w)
    sig = make_signature_from_w(message, salt, w, public_material)

    rejecting_public_material = dict(public_material)
    rejecting_public_material["coset_relation"] = lambda witness, target_h, material: False
    details = verify(message, sig, rejecting_public_material, return_details=True)
    assert_true(not details["valid"], "custom coset relation should be enforced")
    assert_true(not details["checks"]["same_public_coset"], "custom coset relation check did not fail")

    raised = False
    try:
        signature_from_h_and_w(target["h"], w, k=1, public_material=rejecting_public_material)
    except ValueError:
        raised = True
    assert_true(raised, "signature derivation should enforce the public coset relation")


def test_norm_boundary_cases():
    message = "boundary"
    salt = fixed_salt(0x60)
    public_material = make_toy_e8_public_material(k=1)
    target = hash_to_target(message, salt=salt, k=1, public_context=public_material["public_context"], salt_bytes=public_material_salt_bytes(public_material))
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
        test_pinned_basis_and_representative_order,
        test_salt_length_is_fixed,
        test_hawk_parameter_sets,
        test_sampler_rng_hooks,
        test_valid_signature_verifies,
        test_q_public_norm_path,
        test_mismatched_public_material_fails,
        test_material_is_required,
        test_tampered_signature_fails,
        test_negated_signature_fails_by_symbreak,
        test_oversized_witness_fails,
        test_non_e8_witness_fails,
        test_public_coset_relation_is_material_controlled,
        test_norm_boundary_cases,
    ]
    for test in tests:
        test()
        print("PASS", test.__name__)
    print("all E8-HAWK signing-stage tests passed")


run_all_tests()
