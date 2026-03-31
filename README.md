# Numerical Evolution of Closed Strings in pp-Wave Spacetimes

This repository contains a Python-based computational study of closed string dynamics in the background of a gravitational plane wave (pp-wave). The project simulates the interaction between a string and various gravitational pulses to numerically visualize and analyze the **gravitational memory effect**.

## Project Overview

In string theory, pp-waves represent exact solutions for gravitational radiation. When a closed string interacts with such a pulse, its internal degrees of freedom are excited. This simulation solves the equations of motion for the transverse string amplitudes $X^i(\tau, \sigma)$ in the light-cone gauge:

$$\ddot{X}^i + H_{ij}(\tau) X^j = 0$$

where $H_{ij}(\tau)$ defines the profile of the gravitational wave pulse. The "memory effect" is observed as a permanent residual deformation of the string—quantified by its ellipticity—after the pulse has passed.

## Features

* **Numerical Solver:** Integrates time-dependent harmonic oscillator equations for string modes.
* **Pulse Profiles:** Support for multiple pulse shapes:
    * Gaussian
    * $\text{sech}^2$
    * Square wave
* **Analysis Tools:**
    * Calculates time-dependent **ellipticity** to quantify deformation.
    * Visualizes the "stretching and squeezing" of the string in the transverse plane.
* **Visualization:** Generates static plots of string evolution and optional animations of the interaction.

## Technical Requirements

The simulation is written in Python 3.x and requires the following standard scientific libraries:

* `numpy` (Array manipulation and math)
* `matplotlib` (Data visualization and animation)
* `scipy` (Numerical integration via `solve_ivp` or `odeint`)

Install dependencies using:
```bash
pip install numpy matplotlib scipy
