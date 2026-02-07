
# Code for Scrambling for Dollars - Bianchi, Bigio, and Engel 2025

## Overview
This repository contains MATLAB and Julia code for the **Scrambling for Dollars** framework, including empirical filtering and a nonlinear multi-currency model.

## Pipeline

Stage 1: Filtering/main_filter.m        → RW_shock.csv
Stage 2: Filtering/markov_estimation.jl → MS_sigma_us_prob.csv
Stage 3: Filtering/plotting/plot_regimes.m → Figures 2 & 3
Stage 4: Nonlinear Model/main_LFX.m     → Table 4

## Folder Structure

### Filtering/
- main_filter.m — Back out σ shocks from TED spreads
- markov_estimation.jl — Julia Markov estimation
- functions/ — load_data.m, params.m, model_equations.m
- functions/chi/leontief/ — Chi functions (baseline)
- functions/chi/cobb_douglas/ — Chi functions (alternative)
- plotting/ — plot_regimes.m, plot_baseline.m, plot_currencies.m, plot_data.m
- utils/ — exportfig.m, suptitle.m
- data/ — LFX_datainputs.xlsx, LFX_data3.mat, calibration.mat

### Nonlinear Model/
- main_LFX.m — Main model solver
- feqm.m — 7-equation equilibrium system
- feqm_calibrate.m — Calibration targets
- feqm_multicur.m — 2-equation reduced system
- feqm_vec.m — Vectorized equilibrium

## Model Variants
Select chi function specification in main_filter.m:
model_variant = 'leontief';  % or 'cobb_douglas'

## Requirements
- MATLAB (or GNU Octave)
- Julia (for Markov estimation)

## Authors
- Javier Bianchi
- Saki Bigio
- Charles Engel
