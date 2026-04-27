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

from sage.all import ZZ


def sym_break(w):
    # HAWK sym-break analogue: true iff w is nonzero and the first nonzero
    # coordinate is positive.
    for wi in w:
        wi = ZZ(wi)
        if wi > 0:
            return True
        if wi < 0:
            return False
    return False


def canonical_representative(w):
    if sym_break(w):
        return w, False
    minus_w = -w
    if sym_break(minus_w):
        return minus_w, True
    raise ValueError("zero vector has no sym-break representative")
