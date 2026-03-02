# 1D Hunt-Crossley Contact Simulation Design

## Overview

A simple Python simulation to test continuous contact method using the Hunt-Crossley contact model. A ball drops from rest under gravity and bounces on a fixed ground plane.

## Physics Model

### System
- **Body**: Single point mass m in 1D (x direction, positive upward)
- **Ground**: Fixed plane at x=0
- **Penetration**: δ = -x when x < 0

### Hunt-Crossley Contact Force
```
F_contact = k * δ^n * (1 + c * δ_dot)
where δ = -x (penetration), δ_dot = -v (penetration velocity)
```
- k: contact stiffness (N/m^n)
- n: exponent (1.5 for Hertzian-like behavior)
- c: damping coefficient (s/m)

Note: Force is only applied when δ > 0 (penetration exists). The damping term (1 + c * δ_dot) increases force during compression (δ_dot > 0) and decreases it during restitution (δ_dot < 0), providing energy dissipation.

### Equation of Motion
```
m * x'' = -m*g + F_contact   (when x < 0, i.e., penetration)
m * x'' = -m*g               (when x >= 0, free flight)
```

## Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| m | 1.0 kg | Mass |
| g | 9.81 m/s² | Gravity |
| h | 1.0 m | Initial height |
| v₀ | 0 m/s | Initial velocity |
| k | 1e5 N/m^n | Contact stiffness |
| n | 1.5 | Exponent (Hertzian-like) |
| c | 0.1 s/m | Damping coefficient |
| t_end | 3.0 s | Simulation time |

## Implementation

### ODE Function
```python
def equations(t, y):
    x, v = y
    delta = -x  # penetration (positive when x < 0)

    if delta > 0:
        F_contact = k * delta**n * (1 + c * v)
    else:
        F_contact = 0

    a = -g + F_contact / m
    return [v, a]
```

### Solver
- Use `scipy.integrate.solve_ivp` with 'RK45' method
- Set `dense_output=True` for smooth plotting
- Use appropriate tolerances (rtol=1e-6, atol=1e-9)

### Output
- Time array for plotting
- Position and velocity history
- Contact force history (computed from position)

### Visualization
- Main plot: Contact force vs time
- Secondary plot: Position vs time (for reference)

## File Structure

```
continuous_contact/
├── docs/
│   └── plans/
│       └── 2026-03-02-hunt-crossley-contact-design.md
└── drop_test.py
```

## Success Criteria

1. Ball drops from h=1m, bounces multiple times with decreasing height
2. Contact force is zero during free flight, spikes during contact
3. Energy dissipates over time due to damping (visible in decreasing bounce heights)
4. Smooth, physically reasonable plots
