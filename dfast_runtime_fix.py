#!/usr/bin/env python3
"""Run DFAST 1.5.3 with an upstream bug worked around at runtime. Nothing on disk changes.

The bug: dfc/utils/feature_util.py:174 reads

    getattr(hit, "description", hit.id)

Python evaluates a default argument eagerly, so hit.id runs even when hit HAS a
description. HmmHit (dfc/models/hit.py:85) defines accession, name, description and no
id, so this raises AttributeError and aborts the run. It fires only where an rRNA or
tRNA overlaps a CDS named from an HMM hit, which is why 342 of 363 genomes finished and
21 died on exactly this line.

The workaround: give HmmHit a class-level `id` attribute of None, so the eagerly
evaluated fallback has something to read. getattr then still returns `description`,
because HmmHit always sets it in __init__. That line only builds the human-readable
reason string attached to a removal record, so it cannot change which features are
removed or any pseudogene count.

A plain class attribute is used rather than a property: an earlier attempt added a
read-only property to every Hit subclass and broke ProteinHit, which is assigned an
`id` elsewhere in DFAST. A class attribute is shadowed cleanly by any instance
assignment, so nothing else is affected.

Usage, exactly as dfast:
    python dfast_runtime_fix.py -g genome.fna -o outdir --cpu 4 --force
"""
import glob
import os
import runpy
import sys
from pathlib import Path

real = Path(sys.executable).parent / "dfast"
if not real.exists():
    sys.exit(f"dfast not found at {real}")

roots = [r for r in sorted(glob.glob(
            str(Path(sys.executable).parent.parent / "opt" / "dfast-*")))
         if os.path.isdir(os.path.join(r, "dfc"))]
if not roots:
    sys.exit("could not locate the dfast install root containing dfc/")
for r in roots:
    if r not in sys.path:
        sys.path.insert(0, r)

from dfc.models.hit import HmmHit  # noqa: E402

if "id" in vars(HmmHit) or hasattr(HmmHit, "id"):
    print("runtime fix: HmmHit already has id, nothing to do", file=sys.stderr)
else:
    HmmHit.id = None
    print("runtime fix: HmmHit.id = None", file=sys.stderr)

sys.argv = ["dfast"] + sys.argv[1:]
runpy.run_path(str(real), run_name="__main__")
