"""Translate the pinned Mol* real solid harmonic polynomials to NumPy syntax."""

import re
import sys
from pathlib import Path

src = (
    Path(sys.argv[1]) / "src/extensions/alpha-orbitals/spherical-functions.ts"
).read_text()
out = ['"""Real solid harmonics L=0..4, adapted from Mol* (MIT; see NOTICE)."""', ""]
for degree, body in re.findall(r"function L(\d)\([^)]*\)\s*\{([^}]+)\}", src):
    body = (
        body.replace("const ", "")
        .replace(";", "")
        .replace(", yy", "; yy")
        .replace(", zz", "; zz")
        .replace(", yyy", "; yyy")
        .replace(", zzz", "; zzz")
        .replace(", yyyy", "; yyyy")
        .replace(", zzzz", "; zzzz")
    )
    out.append(f"def L{degree}(alpha, x, y, z):\n" + body.rstrip() + "\n")
out.append("FUNCTIONS = (L0, L1, L2, L3, L4)\n")
(
    Path(__file__).resolve().parents[1] / "src/molstar_style_in_pymol/harmonics.py"
).write_text("\n".join(out))
