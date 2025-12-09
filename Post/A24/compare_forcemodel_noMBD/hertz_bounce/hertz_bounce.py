"""
1D ball bounce with Hertzian contact force integrated via RK45.

Outputs (saved alongside this script):
  - bounce_results.csv : columns t, y, v
  - bounce_height.png  : height vs time plot

Run:
    python hertz_bounce.py
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


@dataclass
class BallState:
    """Position/velocity of a point mass moving in 1D (vertical)."""

    y: float  # height [m]
    v: float  # velocity [m/s]
    impacts: List[float] = field(default_factory=list)


def simulate_bounce(
    y0: float = 1.0,
    v0: float = 0.0,
    gravity: float = 9.81,
    k_n: float = 1e5,
    c_n: float = 30.0,
    hertz_exp: float = 1.5,
    max_step: float = 1e-2,
    t_end: float = 2.0,
    rtol: float = 1e-6,
    atol: float = 1e-9,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, BallState]:
    """
    Integrate free-fall with a Hertzian normal contact force at y=0 using RK45.

    Ground force (upwards):
        delta = max(-y, 0)
        F_n = k_n * delta^hertz_exp + c_n * delta_dot (delta_dot = -v)
    No discrete events are used; the stiff nonlinearity captures the bounce.
    """

    def dynamics(_t: float, z: np.ndarray) -> list[float]:
        y, v = z
        delta = max(-y, 0.0)
        if delta > 0.0:
            f_spring = k_n * delta**hertz_exp
            f_damp = c_n * (-v)
            f_contact = f_spring + f_damp
        else:
            f_contact = 0.0
        a = f_contact - gravity  # mass = 1 kg
        return [v, a]

    state = BallState(y=y0, v=v0)

    sol = solve_ivp(
        dynamics,
        (0.0, t_end),
        np.array([y0, v0], dtype=float),
        method="RK45",
        max_step=max_step,
        rtol=rtol,
        atol=atol,
    )

    t = sol.t
    y = sol.y[0]
    v = sol.y[1]

    # Record approximate impact times as local minima when y is near zero and velocity changes sign.
    for i in range(1, len(y) - 1):
        if y[i] <= 0.0 and v[i - 1] < 0.0 <= v[i + 1]:
            state.impacts.append(t[i])

    return t, y, v, state


def main() -> None:
    out_dir = Path(__file__).resolve().parent
    t, y_hist, v_hist, state = simulate_bounce()

    csv_path = out_dir / "bounce_results.csv"
    np.savetxt(csv_path, np.column_stack((t, y_hist, v_hist)), delimiter=",", header="t,y,v", comments="")
    print(f"Saved results to: {csv_path}")

    plt.figure(figsize=(7, 4))
    plt.plot(t, y_hist, label="height y(t)")
    plt.axhline(0.0, color="k", linewidth=1, linestyle="--", label="ground")
    plt.xlabel("time [s]")
    plt.ylabel("height [m]")
    plt.title("1D Ball Bounce (Hertz contact)")
    plt.legend()
    plt.tight_layout()

    png_path = out_dir / "bounce_height.png"
    plt.savefig(png_path, dpi=200)
    print(f"Saved plot to: {png_path}")
    plt.show()


if __name__ == "__main__":
    main()
