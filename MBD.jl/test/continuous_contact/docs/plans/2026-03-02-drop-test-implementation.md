# 1D Hunt-Crossley Drop Test Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Create a Python simulation of a ball dropping and bouncing using the Hunt-Crossley continuous contact model, with contact force visualization.

**Architecture:** Single-file simulation using scipy's solve_ivp for ODE integration. The Hunt-Crossley model provides smooth, continuous contact forces during penetration, making it ideal for testing continuous contact methods in multibody dynamics.

**Tech Stack:** Python 3, numpy, scipy, matplotlib

---

## Task 1: Create Basic Simulation Script

**Files:**
- Create: `drop_test.py`

**Step 1: Create file with imports and parameters**

```python
"""
1D Hunt-Crossley Contact Simulation
A ball drops from rest and bounces on a fixed ground plane.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Physical parameters
m = 1.0          # mass (kg)
g = 9.81         # gravity (m/s^2)
h0 = 1.0         # initial height (m)
v0 = 0.0         # initial velocity (m/s)

# Hunt-Crossley contact parameters
k = 1e5          # contact stiffness (N/m^n)
n = 1.5          # exponent (Hertzian-like)
c = 0.1          # damping coefficient (s/m)

# Simulation parameters
t_end = 3.0      # simulation time (s)
```

**Step 2: Add ODE function**

```python
def hunt_crossley_force(x, v, k, n, c):
    """Calculate Hunt-Crossley contact force."""
    delta = -x  # penetration (positive when x < 0)
    if delta > 0:
        return k * delta**n * (1 + c * v)
    return 0.0


def equations(t, y):
    """ODE system: y = [x, v] where x is position, v is velocity."""
    x, v = y
    F_contact = hunt_crossley_force(x, v, k, n, c)
    a = -g + F_contact / m
    return [v, a]
```

**Step 3: Add solver and solution extraction**

```python
def solve():
    """Solve the ODE and return results."""
    y0 = [h0, v0]  # initial state
    t_span = (0, t_end)
    t_eval = np.linspace(0, t_end, 1000)

    sol = solve_ivp(
        equations,
        t_span,
        y0,
        method='RK45',
        t_eval=t_eval,
        rtol=1e-6,
        atol=1e-9
    )

    return sol.t, sol.y[0], sol.y[1]


def compute_contact_force(x, v):
    """Compute contact force array from position and velocity."""
    return np.array([hunt_crossley_force(xi, vi, k, n, c)
                     for xi, vi in zip(x, v)])
```

**Step 4: Add plotting function**

```python
def plot_results(t, x, v, F):
    """Plot simulation results."""
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    # Position vs time
    axes[0].plot(t, x, 'b-', linewidth=1.5)
    axes[0].axhline(y=0, color='k', linestyle='--', alpha=0.5, label='Ground')
    axes[0].set_ylabel('Position x (m)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # Velocity vs time
    axes[1].plot(t, v, 'g-', linewidth=1.5)
    axes[1].axhline(y=0, color='k', linestyle='--', alpha=0.5)
    axes[1].set_ylabel('Velocity v (m/s)')
    axes[1].grid(True, alpha=0.3)

    # Contact force vs time
    axes[2].plot(t, F, 'r-', linewidth=1.5)
    axes[2].set_ylabel('Contact Force (N)')
    axes[2].set_xlabel('Time (s)')
    axes[2].grid(True, alpha=0.3)

    plt.suptitle('1D Hunt-Crossley Contact Simulation')
    plt.tight_layout()
    plt.show()
```

**Step 5: Add main function**

```python
def main():
    """Main entry point."""
    print("Running 1D Hunt-Crossley contact simulation...")
    print(f"Parameters: m={m} kg, h0={h0} m, k={k:.0e} N/m^{n}, c={c} s/m")

    t, x, v = solve()
    F = compute_contact_force(x, v)

    # Print some statistics
    max_force = np.max(F)
    print(f"Maximum contact force: {max_force:.2f} N")
    print(f"Final position: {x[-1]:.4f} m")
    print(f"Final velocity: {v[-1]:.4f} m/s")

    plot_results(t, x, v, F)


if __name__ == "__main__":
    main()
```

**Step 6: Run the simulation to verify it works**

Run: `cd MBD.jl/test/continuous_contact && python drop_test.py`
Expected: Three plots appear showing position, velocity, and contact force vs time. Ball bounces with decreasing height.

**Step 7: Commit**

```bash
git add MBD.jl/test/continuous_contact/drop_test.py
git add MBD.jl/test/continuous_contact/docs/
git commit -m "feat: add 1D Hunt-Crossley contact simulation for continuous contact testing"
```

---

## Summary

| Task | Description | Files |
|------|-------------|-------|
| 1 | Create simulation script | `drop_test.py` |

Total: 1 file to create, ~100 lines of code.
