# Sampling q-ary lattices

This repository contains samplers some q-ary lattices. We also include code for the samplers in Espitau et. al.'s work following the pseudocode in their paper. 

## Small samplers

Rather than provide the samplers in their most general form, we provide samplers for the root lattices which are used for comparison in the paper. We compare the time it takes to obtain 100,000 samples. 

## Further documentation

We mainly rely on the discrete Gaussian sampler for the integers provided by Sage. Further documentation can be found here: https://doc.sagemath.org/html/en/reference/stats/sage/stats/distributions/discrete_gaussian_integer.html.

## E8 HAWK-style signing stage

The `e8_hawk/` directory contains a modular E8-based prototype of the HAWK signing stage flow. It is not yet a drop in HAWK replacement.

Local preparation already done:

1. `material.sage` defines explicit `secret_material` and `public_material` adapters.
2. `sign.sage` requires external secret material, it does not generate keys.
3. `verify.sage` requires external public material.
4. Hashing is bound to `public_context`, ready for HAWK public key bytes.
5. The verifier norm path supports a public quadratic object `Q` via `w Q w^T`.
6. Tests cover context changes, mismatched material, Q-based norm checks, sym-break, and norm rejection.

Remaining HAWK-specific steps:

1. Map HAWK's real signing key/trapdoor into `secret_material`.
2. Map HAWK's real public key into `public_material`, including `Q` and transcript bytes.
3. Prove that the E8 witness relation being sampled is the relation checked by HAWK verification.
4. Wire `sign.sage`/`verify.sage` into HAWK's encoding, compression, and transcript format.
