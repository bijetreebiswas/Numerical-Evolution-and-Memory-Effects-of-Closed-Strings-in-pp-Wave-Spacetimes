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
  
## Results
<img width="720" height="288" alt="Figure_1" src="https://github.com/user-attachments/assets/143da8ac-bde3-47f0-a7f8-f1f8a7f06703" />
<img width="720" height="284" alt="Figure_2" src="https://github.com/user-attachments/assets/d4a2482d-c396-4390-93dc-f7cc19652912" />
<img width="768" height="377" alt="Figure_3" src="https://github.com/user-attachments/assets/6927d621-ba4f-4007-848b-eb82ab282379" />
<img width="432" height="322" alt="Figure_4" src="https://github.com/user-attachments/assets/e463ab32-9fbf-4529-84b5-9315a10d1759" />

