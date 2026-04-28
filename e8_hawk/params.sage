from sage.all import QQ, ZZ

# Centralized parameters for the E8-based HAWK-style signing stage.
# The named HAWK profiles below use HAWK v1.1 parameter-set dimensions. This
# E8 prototype has direct dimension 8k, while HAWK's signing ambient dimension
# is 2n, so a HAWK degree n maps to k = n/4 E8 blocks.

E8_BLOCK_DIM = ZZ(8)
MODULUS = ZZ(2)

HAWK_PARAMETER_SETS = {
    "HAWK-256": {
        "name": "HAWK-256",
        "degree": ZZ(256),
        "targeted_security": "Challenge",
        "bit_security": ZZ(64),
        "transcript_size_limit_log2": ZZ(32),
        "private_key_bytes": ZZ(96),
        "public_key_bytes": ZZ(450),
        "signature_bytes": ZZ(249),
        "eta": ZZ(2),
        "sigma_sign": QQ(1010) / QQ(1000),
        "sigma_verify": QQ(1042) / QQ(1000),
        "sigma_key_recovery": QQ(1042) / QQ(1000),
        "salt_bits": ZZ(112),
        "keygen_seed_bits": ZZ(128),
        "hashed_pubkey_bits": ZZ(128),
    },
    "HAWK-512": {
        "name": "HAWK-512",
        "degree": ZZ(512),
        "targeted_security": "NIST-I",
        "bit_security": ZZ(128),
        "transcript_size_limit_log2": ZZ(64),
        "private_key_bytes": ZZ(184),
        "public_key_bytes": ZZ(1024),
        "signature_bytes": ZZ(555),
        "eta": ZZ(4),
        "sigma_sign": QQ(1278) / QQ(1000),
        "sigma_verify": QQ(1425) / QQ(1000),
        "sigma_key_recovery": QQ(1425) / QQ(1000),
        "salt_bits": ZZ(192),
        "keygen_seed_bits": ZZ(192),
        "hashed_pubkey_bits": ZZ(256),
    },
    "HAWK-1024": {
        "name": "HAWK-1024",
        "degree": ZZ(1024),
        "targeted_security": "NIST-V",
        "bit_security": ZZ(256),
        "transcript_size_limit_log2": ZZ(64),
        "private_key_bytes": ZZ(360),
        "public_key_bytes": ZZ(2440),
        "signature_bytes": ZZ(1221),
        "eta": ZZ(8),
        "sigma_sign": QQ(1299) / QQ(1000),
        "sigma_verify": QQ(1571) / QQ(1000),
        "sigma_key_recovery": QQ(1974) / QQ(1000),
        "salt_bits": ZZ(320),
        "keygen_seed_bits": ZZ(320),
        "hashed_pubkey_bits": ZZ(512),
    },
}


def hawk_parameter_set(name):
    if name not in HAWK_PARAMETER_SETS:
        raise ValueError("unknown HAWK parameter set")
    params = dict(HAWK_PARAMETER_SETS[name])
    params["ambient_dimension"] = ZZ(2) * params["degree"]
    params["e8_blocks"] = ZZ(params["ambient_dimension"] // E8_BLOCK_DIM)
    params["salt_bytes"] = ZZ(params["salt_bits"] // 8)
    params["keygen_seed_bytes"] = ZZ(params["keygen_seed_bits"] // 8)
    params["hashed_pubkey_bytes"] = ZZ(params["hashed_pubkey_bits"] // 8)
    return params


DEFAULT_HAWK_PARAMETER_SET = "HAWK-512"
HAWK256_BLOCKS = hawk_parameter_set("HAWK-256")["e8_blocks"]
HAWK512_BLOCKS = hawk_parameter_set("HAWK-512")["e8_blocks"]
HAWK1024_BLOCKS = hawk_parameter_set("HAWK-1024")["e8_blocks"]
DEFAULT_BLOCKS = HAWK512_BLOCKS

SALT_BYTES = int(hawk_parameter_set(DEFAULT_HAWK_PARAMETER_SET)["salt_bytes"])
HASH_DOMAIN = b"E8-HAWK-signing-stage-v1"

# Width used by the E8 coset sampler. The sampler draws w from h + 2E8^k
# at width 2 * SIGMA_SIGN, mirroring HAWK's 2*sigma_sign sampling convention.
SIGMA_SIGN = hawk_parameter_set(DEFAULT_HAWK_PARAMETER_SET)["sigma_sign"]
SIGMA_VERIFY = hawk_parameter_set(DEFAULT_HAWK_PARAMETER_SET)["sigma_verify"]

# HAWK's verification bound is 8*n*sigma_verify^2 for dimension 2n.
# For a direct dimension d = 8k, this is 4*d*sigma_verify^2.
NORM_BOUND_FACTOR = QQ(4)
MAX_REJECTION_ATTEMPTS = 1000

# Numerical precision for theta-weight computations in the sampler.
MP_DPS = 80


def e8_dimension(k=DEFAULT_BLOCKS):
    return ZZ(E8_BLOCK_DIM * ZZ(k))


def challenge_bit_count(k=DEFAULT_BLOCKS):
    # One E8/2E8 representative is selected by 8 binary coefficients in the
    # fixed E8 basis; E8^k therefore needs 8k bits.
    return ZZ(E8_BLOCK_DIM * ZZ(k))


def norm_bound_sq(k=DEFAULT_BLOCKS, sigma_verify=SIGMA_VERIFY):
    return QQ(NORM_BOUND_FACTOR) * QQ(e8_dimension(k)) * QQ(sigma_verify) ** 2


def parameter_summary(
    k=None,
    sigma_sign=None,
    sigma_verify=None,
    salt_bytes=None,
    parameter_set=DEFAULT_HAWK_PARAMETER_SET,
):
    params = hawk_parameter_set(parameter_set) if parameter_set is not None else {}
    if k is None:
        k = params.get("e8_blocks", DEFAULT_BLOCKS)
    if sigma_sign is None:
        sigma_sign = params.get("sigma_sign", SIGMA_SIGN)
    if sigma_verify is None:
        sigma_verify = params.get("sigma_verify", SIGMA_VERIFY)
    if salt_bytes is None:
        salt_bytes = params.get("salt_bytes", SALT_BYTES)
    return {
        "k": ZZ(k),
        "dimension": e8_dimension(k),
        "hawk_parameter_set": parameter_set,
        "hawk_degree": params.get("degree", None),
        "hawk_ambient_dimension": params.get("ambient_dimension", None),
        "modulus": MODULUS,
        "sigma_sign": QQ(sigma_sign),
        "sigma_verify": QQ(sigma_verify),
        "norm_bound_sq": norm_bound_sq(k, sigma_verify),
        "salt_bytes": int(ZZ(salt_bytes)),
        "max_rejection_attempts": MAX_REJECTION_ATTEMPTS,
    }
