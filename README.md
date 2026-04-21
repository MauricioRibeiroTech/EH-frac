# Resilience and Optimal Power Flow in Multilayer Networks ⚡🧠
### Dynamics of Fractional Hysteretic Oscillators

[![Julia](https://img.shields.io/badge/Language-Julia-9558B2.svg?style=flat-square&logo=julia)](https://julialang.org/)
[![Python](https://img.shields.io/badge/Language-Python-3776AB.svg?style=flat-square&logo=python)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=flat-square)](https://opensource.org/licenses/MIT)

This repository contains the computational framework developed for the study of **Phase Synchronization** and **Energy Harvesting** in multilayer networks consisting of fractional-order Bouc-Wen oscillators.

The project investigates how complex network topologies (Small-World vs. Random) and the fractional order ($\alpha$) influence the resilience, synchronization stability, and optimal power flow in coupled dynamical systems.

## 🚀 Project Architecture

The framework is organized into two primary computational pipelines:

1.  **Simulation Core (`Julia`):** High-performance implementation using specialized solvers for fractional differential equations and complex network metrics.
    * `2_multicamada.jl`: Core script managing inter-layer piezoelectric coupling ($\sigma_{inter}$) and calculating the Kuramoto Order Parameter ($R$).
2.  **Analysis & Visualization (`Python`):** Scripts optimized for generating publication-quality figures following scientific standards (Nature-style aesthetics).
    * `1_plotting.py`: Processes simulation data to generate comparative dashboards for harvested power and entropy.
    * `2_bianconi.py`: Computes and visualizes Bianconi Structural Entropy across varying network scales ($N$).

## 🧬 Key Features

* **Fractional Calculus:** Integration of hysteretic models with long-term memory effects.
* **Multilayer Analysis:** Mechanical layer (Small-World) coupled with an electrical harvesting layer (Random).
* **Information Theory Metrics:** Calculation of Shannon Entropy and Mutual Information (MI) between non-linear oscillators.
* **Synchronization Dynamics:** Real-time monitoring of phase transitions in large-scale networks.

---
**Developed by Maurício Aparecido Ribeiro** *Researcher in Theoretical Physics and Complex Systems.*
