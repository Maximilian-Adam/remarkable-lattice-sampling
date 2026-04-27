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

from sage.all import ZZ, vector
import hashlib


def normalize_message(message):
    if isinstance(message, bytes):
        return b"bytes" + len(message).to_bytes(8, "big") + message
    if isinstance(message, str):
        data = message.encode("utf-8")
        return b"str" + len(data).to_bytes(8, "big") + data

    bits = []
    try:
        for b in message:
            z = ZZ(b)
            if z not in (0, 1):
                raise ValueError
            bits.append(str(int(z)))
    except Exception as exc:
        raise ValueError("message must be bytes, str, or an iterable of bits") from exc

    bitstring = "".join(bits)
    pad = (-len(bitstring)) % 8
    packed = int(bitstring + ("0" * pad), 2).to_bytes((len(bitstring) + pad) // 8, "big") if bitstring else b""
    return b"bits" + len(bitstring).to_bytes(8, "big") + packed


def normalize_salt(salt=None, salt_bytes=SALT_BYTES):
    if salt is None:
        return os.urandom(salt_bytes)
    if isinstance(salt, bytes):
        return salt
    if isinstance(salt, str):
        s = salt[2:] if salt.startswith("0x") else salt
        if len(s) % 2 != 0:
            s = "0" + s
        try:
            return bytes.fromhex(s)
        except ValueError as exc:
            raise ValueError("salt hex string is invalid") from exc
    raise ValueError("salt must be None, bytes, or a hex string")


def normalize_public_context(public_context=b""):
    if public_context is None:
        return b""
    if isinstance(public_context, bytes):
        return public_context
    if isinstance(public_context, str):
        return public_context.encode("utf-8")
    raise ValueError("public_context must be bytes, str, or None")


def build_hash_transcript(message, salt, k=DEFAULT_BLOCKS, public_context=b"", domain=HASH_DOMAIN):
    msg = normalize_message(message)
    salt_b = normalize_salt(salt)
    ctx = normalize_public_context(public_context)
    return (
        domain
        + int(ZZ(k)).to_bytes(4, "big")
        + len(ctx).to_bytes(4, "big") + ctx
        + len(salt_b).to_bytes(2, "big") + salt_b
        + msg
    )


def shake_256_bits(transcript, bit_count):
    bit_count = ZZ(bit_count)
    if bit_count <= 0:
        raise ValueError("bit_count must be positive")
    digest = hashlib.shake_256(transcript).digest((int(bit_count) + 7) // 8)
    bits = []
    for byte in digest:
        for shift in range(7, -1, -1):
            bits.append(ZZ((byte >> shift) & 1))
            if len(bits) == bit_count:
                return vector(ZZ, bits)
    raise RuntimeError("not enough XOF output")


def hash_to_target(message, salt=None, k=DEFAULT_BLOCKS, public_context=b"", domain=HASH_DOMAIN):
    salt_b = normalize_salt(salt)
    transcript = build_hash_transcript(
        message,
        salt_b,
        k=k,
        public_context=public_context,
        domain=domain,
    )
    bits = shake_256_bits(transcript, challenge_bit_count(k))
    h = coset_rep_from_bits(bits, k=k)
    return {
        "h": h,
        "bits": bits,
        "salt": salt_b,
        "k": ZZ(k),
        "public_context": normalize_public_context(public_context),
        "domain": domain,
        "transcript": transcript,
    }
