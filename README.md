# EH-frac

# Resilience and Optimal Power Flow in Multilayer Networks ⚡🧠
### Dynamics of Fractional Hysteretic Oscillators

[![Julia](https://img.shields.io/badge/Language-Julia-9558B2.svg?style=flat-square&logo=julia)](https://julialang.org/)
[![Python](https://img.shields.io/badge/Language-Python-3776AB.svg?style=flat-square&logo=python)](https://www.python.org/)

This repository provides the computational framework for studying **Phase Synchronization** and **Energy Harvesting** in multilayer networks of fractional-order Bouc-Wen oscillators.

The project investigates how network topology (Small-World vs. Random) and fractional order ($\alpha$) influence the resilience and power transfer efficiency of coupled complex systems.

## 🚀 Project Architecture

The framework is structured into two main pipelines:

1.  **Simulation Core (`Julia`):** High-performance implementation for solving fractional differential equations and network metrics.
    * `2_multicamada.jl`: Main script managing piezoelectric coupling ($\sigma_{inter}$) and calculating the Kuramoto Order Parameter ($R$).
2.  **Analysis & Visualization (`Python`):** Optimized scripts for generating publication-quality figures (Nature Style).
    * `1_plotting.py`: Generates comparative dashboards for power and entropy.
    * `2_bianconi.py`: Visualizes Bianconi Structural Entropy across different network sizes ($N$).

## 📝 Citation

If this framework assists in your research, please cite our work:

> **Ribeiro, M. A.**, et al. (2026). *Resilience and Optimal Power Flow in Multilayer Networks of Fractional Hysteretic Oscillator*. **Nonlinear Dynamics**.
