"""Refresh SNFG reference data from a local Mol* checkout."""

import json
import re
import sys
from pathlib import Path

root = Path(sys.argv[1])
path = Path(__file__).resolve().parents[1] / "src/molstar_style_in_pymol/reference.json"
src = (
    root / "src/mol-model/structure/structure/carbohydrates/constants.ts"
).read_text()
colors = {k: int(v, 16) for k, v in re.findall(r"(\w+): (0x[0-9a-f]+)", src)}
shapes = dict(re.findall(r"\[SaccharideType\.(\w+)\]: SaccharideShape\.(\w+)", src))
components = {}
for abbr, color, kind in re.findall(
    r"abbr: '([^']+)'.*?color: SaccharideColors\.(\w+), type: SaccharideType\.(\w+)",
    src,
):
    components[abbr] = {
        "abbr": abbr,
        "type": kind,
        "shape": shapes[kind],
        "color": colors[color],
    }
result = {}
for key, values in re.findall(
    r"^\s*['\"]?([\w-]+)['\"]?:\s*\[([^]]*)\]", src, re.MULTILINE
):
    if key in components:
        values = re.sub(r"//[^\n]*", "", values)
        for name in re.findall(r"'([^']+)'", values):
            result[name] = components[key]
for key, value in components.items():
    result.setdefault(key.upper(), value)
data = json.loads(path.read_text())
data["saccharides"] = result
path.write_text(json.dumps(data, indent=2) + "\n")
print(len(result), "saccharide residue names")
