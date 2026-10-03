import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from step_graphon import six_block_graphon  # noqa: E402


class SixBlockCandidateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        data = json.loads((ROOT / "data/six_block_candidate.json").read_text())
        cls.graphon = six_block_graphon(data)

    def test_equilibrium(self) -> None:
        self.assertLess(max(map(abs, self.graphon.torques())), 1e-13)

    def test_degree_bottleneck(self) -> None:
        degrees = self.graphon.degrees()
        self.assertAlmostEqual(min(degrees), 0.6915537605671238, places=14)

    def test_rotation_and_two_weak_modes(self) -> None:
        eigenvalues = self.graphon.quotient_hessian_eigenvalues()
        self.assertLess(abs(eigenvalues[0]), 1e-14)
        self.assertGreater(eigenvalues[1], 0.0)
        self.assertGreater(eigenvalues[2], 0.0)
        self.assertLess(eigenvalues[2], 2e-12)

    def test_transverse_step_modes_are_strict(self) -> None:
        self.assertGreater(min(self.graphon.transverse_hessian_values()), 0.24)


if __name__ == "__main__":
    unittest.main()

