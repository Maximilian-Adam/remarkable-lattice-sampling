from sage.all import QQ, ZZ

# Centralized parameters for the E8-based HAWK-style signing stage.
# The tests use small k values for speed, experiments can set k=64 to match
# a 512-dimensional ambient space.

E8_BLOCK_DIM = ZZ(8)
DEFAULT_BLOCKS = ZZ(1)
HAWK512_BLOCKS = ZZ(64)

MODULUS = ZZ(2)
SALT_BYTES = 16
HASH_DOMAIN = b"E8-HAWK-signing-stage-v1"

# Width used by the E8 coset sampler. The sampler draws w from h + 2E8^k
# at width 2 * SIGMA_SIGN, mirroring HAWK's 2*sigma_sign sampling convention.
SIGMA_SIGN = QQ(11) / QQ(5)      # 2.2
SIGMA_VERIFY = QQ(11) / QQ(5)    # 2.2

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


def parameter_summary(k=DEFAULT_BLOCKS, sigma_sign=SIGMA_SIGN, sigma_verify=SIGMA_VERIFY):
    return {
        "k": ZZ(k),
        "dimension": e8_dimension(k),
        "modulus": MODULUS,
        "sigma_sign": QQ(sigma_sign),
        "sigma_verify": QQ(sigma_verify),
        "norm_bound_sq": norm_bound_sq(k, sigma_verify),
        "salt_bytes": SALT_BYTES,
        "max_rejection_attempts": MAX_REJECTION_ATTEMPTS,
    }
