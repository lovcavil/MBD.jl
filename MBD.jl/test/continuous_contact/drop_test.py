"""
1D Contact Simulation - Hunt-Crossley vs Flores Method
Compare different contact force models for continuous contact behavior.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Physical parameters
m = 1.0          # mass (kg)
g = 9.81         # gravity (m/s^2)
h0 = 0.0         # initial height (m)
v0 = -1.0        # initial velocity (m/s)

# Contact parameters
k = 1e5          # contact stiffness (N/m^n)
n = 1.5          # exponent (Hertzian-like)

# Hunt-Crossley damping coefficients to compare
c_values = [0.0, 0.2, 0.5, 1.0]  # damping coefficient (s/m)

# Flores coefficient of restitution to compare
Cr_values = [1.0, 0.5, 0.4, 0.3]  # coefficient of restitution (1.0 = elastic)

# Simulation parameters
t_end = 0.1      # simulation time (s)


def hunt_crossley_force(x, v, k, n, c):
    """Calculate Hunt-Crossley contact force.

    F = k * δ^n * (1 + c * δ_dot)
    where δ = -x (penetration), δ_dot = -v (penetration velocity)
    """
    delta = -x
    delta_dot = -v
    if delta > 0:
        return k * delta**n * (1 + c * delta_dot)
    return 0.0


def flores_force(x, v, k, n, Cr, init_vel):
    """Calculate Flores contact force.

    F = k * δ^n * (1 + eta * v)
    where eta = 8 * (1 - Cr) / (5 * Cr * |init_vel|)

    Note: Force is clamped to be non-negative (tensile force not allowed)
    """
    delta = -x
    if delta > 0:
        # Flores damping coefficient
        eta = 8 * (1.0 - Cr) / (5.0 * Cr * init_vel)
        force = k * delta**n * (1 + eta * v)
        return max(0.0, force)  # No tensile (pulling) force allowed
    return 0.0


def solve_hunt_crossley(c):
    """Solve using Hunt-Crossley method."""
    def equations(t, y):
        x, v = y
        F_contact = hunt_crossley_force(x, v, k, n, c)
        a = -g + F_contact / m
        return [v, a]

    sol = solve_ivp(
        equations, (0, t_end), [h0, v0],
        method='RK45', t_eval=np.linspace(0, t_end, 1000),
        rtol=1e-6, atol=1e-9
    )

    x, v = sol.y[0], sol.y[1]
    F = np.array([hunt_crossley_force(xi, vi, k, n, c) for xi, vi in zip(x, v)])
    return sol.t, x, v, F


def solve_flores(Cr):
    """Solve using Flores method."""
    def equations(t, y):
        x, v = y
        F_contact = flores_force(x, v, k, n, Cr, v0)
        a = -g + F_contact / m
        return [v, a]

    sol = solve_ivp(
        equations, (0, t_end), [h0, v0],
        method='RK45', t_eval=np.linspace(0, t_end, 1000),
        rtol=1e-6, atol=1e-9
    )

    x, v = sol.y[0], sol.y[1]
    F = np.array([flores_force(xi, vi, k, n, Cr, v0) for xi, vi in zip(x, v)])
    return sol.t, x, v, F


def plot_comparison(hc_results, flores_results):
    """Plot comparison of Hunt-Crossley and Flores methods."""

    # Figure 1: Time history comparison - Hunt-Crossley
    fig1, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    colors_hc = plt.cm.Blues(np.linspace(0.4, 0.9, len(c_values)))

    for (c, t, x, v, F), color in zip(hc_results, colors_hc):
        label = f'c = {c} s/m'
        axes[0].plot(t, x, color=color, linewidth=1.5, label=label)
        axes[1].plot(t, v, color=color, linewidth=1.5, label=label)
        axes[2].plot(t, F, color=color, linewidth=1.5, label=label)

    axes[0].axhline(y=0, color='k', linestyle='--', alpha=0.5)
    axes[0].set_ylabel('Position x (m)')
    axes[0].legend(loc='upper right')
    axes[0].grid(True, alpha=0.3)
    axes[1].axhline(y=0, color='k', linestyle='--', alpha=0.5)
    axes[1].set_ylabel('Velocity v (m/s)')
    axes[1].grid(True, alpha=0.3)
    axes[2].set_ylabel('Contact Force (N)')
    axes[2].set_xlabel('Time (s)')
    axes[2].grid(True, alpha=0.3)
    fig1.suptitle('Hunt-Crossley Method - Damping Comparison')
    fig1.tight_layout()
    fig1.savefig('hunt_crossley_time.png', dpi=150)
    print("Saved: hunt_crossley_time.png")

    # Figure 2: Time history comparison - Flores
    fig2, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    colors_fl = plt.cm.Oranges(np.linspace(0.4, 0.9, len(Cr_values)))

    for (Cr, t, x, v, F), color in zip(flores_results, colors_fl):
        label = f'Cr = {Cr}'
        axes[0].plot(t, x, color=color, linewidth=1.5, label=label)
        axes[1].plot(t, v, color=color, linewidth=1.5, label=label)
        axes[2].plot(t, F, color=color, linewidth=1.5, label=label)

    axes[0].axhline(y=0, color='k', linestyle='--', alpha=0.5)
    axes[0].set_ylabel('Position x (m)')
    axes[0].legend(loc='upper right')
    axes[0].grid(True, alpha=0.3)
    axes[1].axhline(y=0, color='k', linestyle='--', alpha=0.5)
    axes[1].set_ylabel('Velocity v (m/s)')
    axes[1].grid(True, alpha=0.3)
    axes[2].set_ylabel('Contact Force (N)')
    axes[2].set_xlabel('Time (s)')
    axes[2].grid(True, alpha=0.3)
    fig2.suptitle('Flores Method - Coefficient of Restitution Comparison')
    fig2.tight_layout()
    fig2.savefig('flores_time.png', dpi=150)
    print("Saved: flores_time.png")

    # Figure 3: Hysteresis comparison - both methods
    fig3, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Hunt-Crossley hysteresis
    for (c, t, x, v, F), color in zip(hc_results, colors_hc):
        delta = -x
        mask = delta > 0
        if np.any(mask):
            axes[0].plot(delta[mask]*1000, F[mask], color=color, linewidth=1.5, label=f'c={c}')
    axes[0].set_xlabel('Indentation (mm)')
    axes[0].set_ylabel('Contact Force (N)')
    axes[0].set_title('Hunt-Crossley: Force vs Indentation')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # Flores hysteresis
    for (Cr, t, x, v, F), color in zip(flores_results, colors_fl):
        delta = -x
        mask = delta > 0
        if np.any(mask):
            axes[1].plot(delta[mask]*1000, F[mask], color=color, linewidth=1.5, label=f'Cr={Cr}')
    axes[1].set_xlabel('Indentation (mm)')
    axes[1].set_ylabel('Contact Force (N)')
    axes[1].set_title('Flores: Force vs Indentation')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    fig3.tight_layout()
    fig3.savefig('hysteresis_comparison.png', dpi=150)
    print("Saved: hysteresis_comparison.png")

    # Figure 4: Direct method comparison (select one from each)
    fig4, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Select c=0.5 and Cr=0.5 for comparison
    hc_selected = hc_results[2]  # c=0.5
    fl_selected = flores_results[2]  # Cr=0.5

    # Position
    axes[0, 0].plot(hc_selected[1], hc_selected[2], 'b-', linewidth=1.5, label=f'Hunt-Crossley c=0.5')
    axes[0, 0].plot(fl_selected[1], fl_selected[2], 'r--', linewidth=1.5, label=f'Flores Cr=0.5')
    axes[0, 0].axhline(y=0, color='k', linestyle=':', alpha=0.5)
    axes[0, 0].set_xlabel('Time (s)')
    axes[0, 0].set_ylabel('Position (m)')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].set_title('Position Comparison')

    # Velocity
    axes[0, 1].plot(hc_selected[1], hc_selected[3], 'b-', linewidth=1.5, label='Hunt-Crossley c=0.5')
    axes[0, 1].plot(fl_selected[1], fl_selected[3], 'r--', linewidth=1.5, label='Flores Cr=0.5')
    axes[0, 1].axhline(y=0, color='k', linestyle=':', alpha=0.5)
    axes[0, 1].set_xlabel('Time (s)')
    axes[0, 1].set_ylabel('Velocity (m/s)')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].set_title('Velocity Comparison')

    # Force
    axes[1, 0].plot(hc_selected[1], hc_selected[4], 'b-', linewidth=1.5, label='Hunt-Crossley c=0.5')
    axes[1, 0].plot(fl_selected[1], fl_selected[4], 'r--', linewidth=1.5, label='Flores Cr=0.5')
    axes[1, 0].set_xlabel('Time (s)')
    axes[1, 0].set_ylabel('Force (N)')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].set_title('Contact Force Comparison')

    # Hysteresis overlay
    delta_hc = -hc_selected[2]
    delta_fl = -fl_selected[2]
    mask_hc = delta_hc > 0
    mask_fl = delta_fl > 0
    axes[1, 1].plot(delta_hc[mask_hc]*1000, hc_selected[4][mask_hc], 'b-', linewidth=1.5, label='Hunt-Crossley c=0.5')
    axes[1, 1].plot(delta_fl[mask_fl]*1000, fl_selected[4][mask_fl], 'r--', linewidth=1.5, label='Flores Cr=0.5')
    axes[1, 1].set_xlabel('Indentation (mm)')
    axes[1, 1].set_ylabel('Force (N)')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].set_title('Hysteresis Comparison')

    fig4.suptitle('Direct Method Comparison: Hunt-Crossley vs Flores')
    fig4.tight_layout()
    fig4.savefig('method_comparison.png', dpi=150)
    print("Saved: method_comparison.png")

    plt.show()


def main():
    """Main entry point."""
    print("=" * 70)
    print("1D Contact Simulation - Hunt-Crossley vs Flores Method")
    print("=" * 70)
    print(f"Parameters: m={m} kg, h0={h0} m, v0={v0} m/s")
    print(f"            k={k:.0e} N/m^{n}")
    print()

    # Hunt-Crossley results
    print("-" * 70)
    print("Hunt-Crossley Method")
    print(f"Damping coefficients: {c_values}")
    hc_results = []
    for c in c_values:
        t, x, v, F = solve_hunt_crossley(c)
        hc_results.append((c, t, x, v, F))
        max_force = np.max(F)
        max_indent = np.max(-x) * 1000
        print(f"  c = {c:.1f} s/m: Max Force = {max_force:8.2f} N, Max Indent = {max_indent:6.3f} mm")

    # Flores results
    print()
    print("-" * 70)
    print("Flores Method")
    print(f"Coefficient of restitution: {Cr_values}")
    flores_results = []
    for Cr in Cr_values:
        t, x, v, F = solve_flores(Cr)
        flores_results.append((Cr, t, x, v, F))
        max_force = np.max(F)
        max_indent = np.max(-x) * 1000
        eta = 8 * (1.0 - Cr) / (5.0 * Cr * v0)
        print(f"  Cr = {Cr:.1f}: eta = {eta:.3f}, Max Force = {max_force:8.2f} N, Max Indent = {max_indent:6.3f} mm")

    print()
    plot_comparison(hc_results, flores_results)


if __name__ == "__main__":
    main()
