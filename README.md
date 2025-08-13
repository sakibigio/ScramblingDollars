# Code for Scrambling  for Dollars - Bianchi, Bigio, and Engel 2025

## Overview
This repository contains MATLAB code for a nonlinear multi-currency model in the **Scrambling for Dollars** framework. The model captures interactions between two currencies (e.g., Dollars and Euros) and includes conditions for reserve demand, bond and deposit, pricing and the determination of the exchange rate.

## Key Files
- **Nonlinear Model/feqm.m** – Baseline seven-equation system describing equilibrium in the two-currency model, linking domestic and foreign interest rates, currency shares, and policy parameters.
- **Nonlinear Model/feqm_calibrate.m** – Extends the baseline system to include calibration targets such as interest rate differentials, exchange rates, and money share dynamics.
- **Nonlinear Model/feqm_multicur.m** – Reduced two-equation system useful for exploring comparative statics of currency demand under different policy environments.
- **Nonlinear Model/feqm_vec.m** – Vectorized version of the equilibrium system to evaluate multiple states simultaneously.

The `Nonlinear Model` directory also contains `.mat` data files and helper routines used in solving and calibrating the model.

## Usage
1. Open MATLAB (or GNU Octave).
2. Navigate to the `Nonlinear Model` directory.
3. Define or supply the helper functions `Echi_d` and `Echi_m` that compute expected marginal utilities of cash and money holdings.
4. Use a nonlinear solver (e.g., `fsolve`) with `feqm.m` or `feqm_calibrate.m` to solve for equilibrium variables.
5. Adjust parameters or targets to reproduce different scenarios, or employ the vectorized version to compute equilibria across states.

## Requirements
- MATLAB or compatible environment such as GNU Octave.

## License
No explicit license is provided. Please contact the original authors for usage permissions.

