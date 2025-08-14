# Design and Simulate the Aerodynamics of Propellers

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository contains my study notes, code, and simulation tools from the Udemy course **“Design and Simulate the Aerodynamics of Propellers”**. The project covers aerodynamic analysis methods, including:

- **Momentum Theory (MT)**
- **Blade Element Momentum Theory (BEMT)**
- **Comparisons with MCEVS (Multi-Configurational E-Vtol Sizing)**

---

## 📂 Repository Structure

- `BEMT.py` &mdash; Python implementation of Blade Element Momentum Theory for propeller performance prediction, with and without MCEVS.
- `MT.py` &mdash; Python implementation of classic Momentum Theory and MCEVS-based calculations for hover climb/descent power.
- `theory.pdf` &mdash; Course notes and theoretical background on propeller aerodynamics, including derivations and formulae.
- `aerodata.xlsx` &mdash; Contains propeller geometry (section radius, pitch, chord, etc.) and airfoil performance data (CL, CD vs. alpha) for use in simulations and interpolation.

---

## 🚀 Features

- **Simulation of propeller performance** in different flight regimes (takeoff, climb, cruise).
- **Comparison of classical methods** (MT, BEMT) with advanced models (MCEVS) that include tip/root loss corrections.
- **Visualization** of thrust, torque, power, efficiency, and dimensionless coefficients.
- **Reusable code** for academic or personal study.
