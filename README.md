# Design and Simulate the Aerodynamics of Propellers

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository contains my study notes and Python implementations of the MATLAB code from the Udemy course **“Design and Simulate the Aerodynamics of Propellers in MATLAB”**.

[View my Udemy certificate](https://www.udemy.com/certificate/UC-4dfa6f8a-e50c-4236-892a-2cda2a33b01a/)

The project covers aerodynamic analysis methods, including:

- **Momentum Theory (MT)**
- **Blade Element Momentum Theory (BEMT)**
- **Comparisons with MCEVS (Multi-Configurational E-Vtol Sizing)**

---

## 📂 Repository Structure

- `BEMT.py` &mdash; Python implementation of Blade Element Momentum Theory for propeller performance prediction.
- `MT.py` &mdash; Python implementation of classic Momentum Theory calculations for hover climb/descent power.
- `theory.pdf` &mdash; Course notes and theoretical background on propeller aerodynamics, including derivations and formulae.
- `aerodata.xlsx` &mdash; Contains propeller geometry (section radius, pitch, chord, etc.) and airfoil performance data (CL, CD vs. alpha) for use in simulations and interpolation.

---

## 🚀 Features

- **Simulation of propeller performance** in different flight regimes (takeoff, climb, cruise).
- **Comparison of methods in the course** (MT, BEMT) with my models (MCEVS) that include tip/root loss corrections.
- **Visualization** of thrust, torque, power, efficiency, and dimensionless coefficients.
- **Reusable code** for academic or personal study.

![BEMT Results](https://github.com/alfiyandyhr/propeller_aerodynamics/raw/main/BEMT_results.png)
