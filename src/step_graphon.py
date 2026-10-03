"""Finite-step Kuramoto graphon calculations using only the Python standard library."""

from __future__ import annotations

from dataclasses import dataclass
import math


def symmetric_eigenvalues(matrix: list[list[float]], tolerance: float = 1e-15) -> list[float]:
    """Jacobi eigenvalue algorithm for a small real symmetric matrix."""
    a = [row[:] for row in matrix]
    n = len(a)
    for _ in range(100 * n * n):
        p, q = max(
            ((i, j) for i in range(n) for j in range(i + 1, n)),
            key=lambda ij: abs(a[ij[0]][ij[1]]),
        )
        if abs(a[p][q]) < tolerance:
            break
        tau = (a[q][q] - a[p][p]) / (2.0 * a[p][q])
        t = math.copysign(1.0, tau) / (abs(tau) + math.sqrt(1.0 + tau * tau)) if tau else 1.0
        cosine = 1.0 / math.sqrt(1.0 + t * t)
        sine = t * cosine
        app, aqq, apq = a[p][p], a[q][q], a[p][q]
        a[p][p] = cosine * cosine * app - 2.0 * sine * cosine * apq + sine * sine * aqq
        a[q][q] = sine * sine * app + 2.0 * sine * cosine * apq + cosine * cosine * aqq
        a[p][q] = a[q][p] = 0.0
        for k in range(n):
            if k in (p, q):
                continue
            akp, akq = a[k][p], a[k][q]
            a[k][p] = a[p][k] = cosine * akp - sine * akq
            a[k][q] = a[q][k] = sine * akp + cosine * akq
    return sorted(a[i][i] for i in range(n))


@dataclass(frozen=True)
class StepGraphon:
    masses: tuple[float, ...]
    weights: tuple[tuple[float, ...], ...]
    phases: tuple[float, ...]

    def validate(self, tolerance: float = 1e-12) -> None:
        n = len(self.masses)
        if len(self.phases) != n or len(self.weights) != n:
            raise ValueError("masses, phases, and weight matrix dimensions differ")
        if any(len(row) != n for row in self.weights):
            raise ValueError("weight matrix is not square")
        if any(p <= 0.0 for p in self.masses):
            raise ValueError("cell masses must be positive")
        if abs(sum(self.masses) - 1.0) > tolerance:
            raise ValueError("cell masses do not sum to one")
        for i in range(n):
            for j in range(n):
                if not 0.0 <= self.weights[i][j] <= 1.0:
                    raise ValueError("graphon weights must lie in [0,1]")
                if abs(self.weights[i][j] - self.weights[j][i]) > tolerance:
                    raise ValueError("graphon weight matrix is not symmetric")

    def torques(self) -> list[float]:
        return [
            sum(
                self.masses[j]
                * self.weights[i][j]
                * math.sin(self.phases[j] - self.phases[i])
                for j in range(len(self.masses))
            )
            for i in range(len(self.masses))
        ]

    def degrees(self) -> list[float]:
        return [
            sum(self.masses[j] * self.weights[i][j] for j in range(len(self.masses)))
            for i in range(len(self.masses))
        ]

    def transverse_hessian_values(self) -> list[float]:
        """The multiplication coefficients q_i on cell-zero-mean modes."""
        return [
            sum(
                self.masses[j]
                * self.weights[i][j]
                * math.cos(self.phases[i] - self.phases[j])
                for j in range(len(self.masses))
            )
            for i in range(len(self.masses))
        ]

    def quotient_hessian(self) -> list[list[float]]:
        n = len(self.masses)
        result = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                cosine_weight = self.weights[i][j] * math.cos(self.phases[i] - self.phases[j])
                result[i][i] += self.masses[j] * cosine_weight
                result[i][j] = -math.sqrt(self.masses[i] * self.masses[j]) * cosine_weight
        return result

    def quotient_hessian_eigenvalues(self) -> list[float]:
        return symmetric_eigenvalues(self.quotient_hessian())


def six_block_graphon(data: dict[str, object]) -> StepGraphon:
    x, y = float(data["x"]), float(data["y"])
    weights = (
        (1.0, x, y, 0.0, 1.0, 1.0),
        (x, 1.0, 0.0, y, 1.0, 1.0),
        (y, 0.0, 1.0, 1.0, 0.0, 1.0),
        (0.0, y, 1.0, 1.0, 1.0, 0.0),
        (1.0, 1.0, 0.0, 1.0, 1.0, 0.0),
        (1.0, 1.0, 1.0, 0.0, 0.0, 1.0),
    )
    graphon = StepGraphon(
        tuple(map(float, data["masses"])),
        weights,
        tuple(map(float, data["phases"])),
    )
    graphon.validate()
    return graphon

