#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from step_graphon import six_block_graphon  # noqa: E402


def main() -> None:
    data = json.loads((ROOT / "data/six_block_candidate.json").read_text())
    graphon = six_block_graphon(data)
    eigenvalues = graphon.quotient_hessian_eigenvalues()
    result = {
        "maximum_absolute_torque": max(map(abs, graphon.torques())),
        "degrees": graphon.degrees(),
        "minimum_degree": min(graphon.degrees()),
        "quotient_hessian_eigenvalues": eigenvalues,
        "two_weak_nonrotation_modes": eigenvalues[1:3],
        "transverse_hessian_values": graphon.transverse_hessian_values(),
        "strictly_positive_transverse_values": min(graphon.transverse_hessian_values()) > 0.0,
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

