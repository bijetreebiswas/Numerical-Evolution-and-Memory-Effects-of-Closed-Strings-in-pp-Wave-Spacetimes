# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 21:02:46 2026

@author: bijet
"""

"""
Numerical simulation of closed string evolution in a pp-wave spacetime
with an arbitrary gravitational wave pulse shape.


Date: 2026-03-30
"""

import matplotlib
matplotlib.use('TkAgg')   # or 'Qt5Agg' if you have PyQt5 installed

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation
import shutil

# ============================================================================
# Pulse shape definitions
# ============================================================================

def gaussian_pulse(tau, A, mu, sigma):
    """Gaussian pulse: A * exp(-(tau-mu)^2/(2*sigma^2))"""
    return A * np.exp(-0.5 * ((tau - mu) / sigma) ** 2)

def sech2_pulse(tau, A, alpha):
    """Sech-squared pulse: A * sech^2(alpha * tau)"""
    return A * (1.0 / np.cosh(alpha * tau)) ** 2

def square_pulse(tau, A, t_start, t_end):
    """Square pulse: A for t_start <= tau <= t_end, else 0."""
    return A if t_start <= tau <= t_end else 0.0

# ============================================================================
# ODE system: convert second-order to first-order
# ============================================================================

def string_ode_system(tau, y, k1, p, W_func, W_args):
    """
    ODE system for [X, dX/dtau, Y, dY/dtau].
    """
    X, Xdot, Y, Ydot = y
    W_val = W_func(tau, *W_args)  # evaluate pulse at current tau

    # Equations:
    # Xddot = - (k1^2 - p^2 * W) * X
    # Yddot = - (k1^2 + p^2 * W) * Y
    Xddot = - (k1**2 - p**2 * W_val) * X
    Yddot = - (k1**2 + p**2 * W_val) * Y

    return [Xdot, Xddot, Ydot, Yddot]

# ============================================================================
# Solve for X(tau) and Y(tau) over a tau range
# ============================================================================

def solve_string_amplitudes(tau_span, tau_eval, initial_conditions, k1, p,
                            W_func, W_args):
    """
    Solve the ODEs for X(tau) and Y(tau) on the interval [tau_span[0], tau_span[1]].
    """
    sol = solve_ivp(
        fun=lambda tau, y: string_ode_system(tau, y, k1, p, W_func, W_args),
        t_span=tau_span,
        y0=initial_conditions,
        t_eval=tau_eval,
        method='RK45',
        rtol=1e-9,
        atol=1e-12
    )
    if not sol.success:
        raise RuntimeError("Integration failed: " + sol.message)
    return sol.y[0], sol.y[2]   # X and Y

# ============================================================================
# Compute string shape at a given tau
# ============================================================================

def string_shape(X, Y, sigma_grid, k1):
    """
    Generate x,y coordinates of the string at a given tau.
    """
    x_vals = X * np.cos(k1 * sigma_grid)
    y_vals = Y * np.sin(k1 * sigma_grid)
    return x_vals, y_vals

# ============================================================================
# Check if ffmpeg is available
# ============================================================================

def is_ffmpeg_available():
    """Return True if ffmpeg executable is found in PATH."""
    return shutil.which('ffmpeg') is not None

# ============================================================================
# Main simulation and plotting
# ============================================================================

def main():
    # ---------- Parameters ----------
    R = 1.0           # initial amplitude (radius of the circle at tau->-inf)
    k1 = 1.0          # mode number
    p = 1.0           # from u = p*tau (can absorb into W)

    # Pulse settings (choose one):
    pulse_type = 'gaussian'    # 'gaussian', 'sech2', or 'square'

    # Set pulse parameters based on type
    if pulse_type == 'gaussian':
        A = 2.0                # amplitude
        mu = 0.0               # center
        sigma = 1.0            # width
        W_func = gaussian_pulse
        W_args = (A, mu, sigma)

    elif pulse_type == 'sech2':
        A = 2.0                # amplitude
        alpha = 1.0            # width parameter
        W_func = sech2_pulse
        W_args = (A, alpha)

    elif pulse_type == 'square':
        A = 2.0                # height
        t_start = -1.0
        t_end = 1.0
        W_func = square_pulse
        W_args = (A, t_start, t_end)
    else:
        raise ValueError("Unknown pulse_type")

    # Time range for integration
    tau_initial = -10.0
    tau_final = 10.0
    tau_span = (tau_initial, tau_final)

    # Number of points for output
    n_tau = 500
    tau_eval = np.linspace(tau_initial, tau_final, n_tau)

    # Initial conditions at tau = tau_initial (circular string in the past)
    X0 = R * np.cos(k1 * tau_initial)
    Xdot0 = -R * k1 * np.sin(k1 * tau_initial)
    Y0 = R * np.cos(k1 * tau_initial)
    Ydot0 = -R * k1 * np.sin(k1 * tau_initial)
    initial_conditions = [X0, Xdot0, Y0, Ydot0]

    # Solve ODEs
    print("Solving ODEs...")
    X_amp, Y_amp = solve_string_amplitudes(tau_span, tau_eval,
                                           initial_conditions,
                                           k1, p, W_func, W_args)
    print("Done.")

    # ---------- Plot pulse shape ----------
    plt.figure(figsize=(10, 4))
    W_vals = np.array([W_func(t, *W_args) for t in tau_eval])
    plt.plot(tau_eval, W_vals, 'b-', linewidth=2)
    plt.xlabel(r'$\tau$', fontsize=12)
    plt.ylabel(r'$W(\tau)$', fontsize=12)
    plt.title(f'Pulse shape: {pulse_type}', fontsize=14)
    plt.grid(True)
    plt.show()

    # ---------- Plot amplitudes X(tau) and Y(tau) ----------
    plt.figure(figsize=(10, 6))
    plt.plot(tau_eval, X_amp, 'r-', label='X(tau)')
    plt.plot(tau_eval, Y_amp, 'b-', label='Y(tau)')
    plt.xlabel(r'$\tau$', fontsize=12)
    plt.ylabel('Amplitude', fontsize=12)
    plt.legend(fontsize=12)
    plt.title('String amplitudes before and after the pulse', fontsize=14)
    plt.grid(True)
    plt.show()

    # ---------- Plot string shape at selected times ----------
    sigma_grid = np.linspace(0, 2*np.pi, 200)
    tau_before = tau_initial + 0.5
    tau_during = 0.0
    tau_after = tau_final - 0.5

    idx_before = np.argmin(np.abs(tau_eval - tau_before))
    idx_during = np.argmin(np.abs(tau_eval - tau_during))
    idx_after = np.argmin(np.abs(tau_eval - tau_after))

    X_before, Y_before = X_amp[idx_before], Y_amp[idx_before]
    X_during, Y_during = X_amp[idx_during], Y_amp[idx_during]
    X_after, Y_after = X_amp[idx_after], Y_amp[idx_after]

    x_before, y_before = string_shape(X_before, Y_before, sigma_grid, k1)
    x_during, y_during = string_shape(X_during, Y_during, sigma_grid, k1)
    x_after, y_after = string_shape(X_after, Y_after, sigma_grid, k1)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    axes[0].plot(x_before, y_before, 'k-', linewidth=2)
    axes[0].set_aspect('equal')
    axes[0].set_title(f'Before pulse ($\\tau = {tau_before:.1f}$)', fontsize=12)
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    axes[0].grid(True)

    axes[1].plot(x_during, y_during, 'k-', linewidth=2)
    axes[1].set_aspect('equal')
    axes[1].set_title(f'During pulse ($\\tau = {tau_during:.1f}$)', fontsize=12)
    axes[1].set_xlabel('x')
    axes[1].grid(True)

    axes[2].plot(x_after, y_after, 'k-', linewidth=2)
    axes[2].set_aspect('equal')
    axes[2].set_title(f'After pulse ($\\tau = {tau_after:.1f}$)', fontsize=12)
    axes[2].set_xlabel('x')
    axes[2].grid(True)

    plt.tight_layout()
    plt.show()

    # ---------- Ellipticity after pulse ----------
    W_all = np.array([W_func(t, *W_args) for t in tau_eval])
    peak_W = np.max(np.abs(W_all))
    threshold = 0.01 * peak_W if peak_W > 0 else 0
    after_idx = (np.abs(W_all) < threshold) & (tau_eval > tau_initial)
    if not np.any(after_idx):
        after_idx = tau_eval > (tau_initial + 0.8 * (tau_final - tau_initial))

    X_after_vals = X_amp[after_idx]
    Y_after_vals = Y_amp[after_idx]
    max_X = np.max(np.abs(X_after_vals))
    max_Y = np.max(np.abs(Y_after_vals))
    ellipticity = max_X / max_Y if max_Y != 0 else np.inf
    print(f"\nEllipticity after pulse (max|X|/max|Y|) = {ellipticity:.4f}")

    # ---------- Animation ----------
    make_animation = input("\nCreate animation? (y/n): ").strip().lower()
    if make_animation == 'y':
        print("Preparing animation...")
        fig_anim, ax_anim = plt.subplots(figsize=(6,6))
        ax_anim.set_aspect('equal')
        ax_anim.set_xlim(-1.5, 1.5)
        ax_anim.set_ylim(-1.5, 1.5)
        ax_anim.grid(True)
        line, = ax_anim.plot([], [], 'b-', linewidth=2)
        title = ax_anim.set_title('')

        def init():
            line.set_data([], [])
            title.set_text('')
            return line, title

        def update(frame):
            tau = tau_eval[frame]
            X = X_amp[frame]
            Y = Y_amp[frame]
            x, y = string_shape(X, Y, sigma_grid, k1)
            line.set_data(x, y)
            title.set_text(f'$\\tau = {tau:.2f}$')
            return line, title

        anim = FuncAnimation(fig_anim, update, frames=len(tau_eval),
                             init_func=init, blit=False, interval=50)

        # Try to display; if it fails, offer save
        try:
            plt.show(block=True)
        except Exception as e:
            print(f"Interactive display failed: {e}")
            save = input("Save as GIF instead? (y/n): ").strip().lower()
            if save == 'y':
                anim.save('string_evolution.gif', writer='pillow')
                print("Saved as string_evolution.gif")
            else:
                print("Animation not saved.")

    print("\nSimulation completed.")

if __name__ == "__main__":
    main()