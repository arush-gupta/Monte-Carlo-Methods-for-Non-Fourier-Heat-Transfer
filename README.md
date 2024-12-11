# Monte Carlo Simulations for Non-Fourier Heat Transfer Using Boltzmann Transport Equation

This repository contains the code and report for a research project focused on modeling heat transfer at the nanoscale, where Fourier's law of heat conduction no longer holds. The project utilizes Monte Carlo simulations to model phonon transport at the nanoscale, an important aspect of understanding thermal conductivity in nanomaterials.

## Overview

At the nanoscale, phonon transport deviates from the classical diffusive behavior described by Fourier's law. Instead, phonons exhibit ballistic or quasi-ballistic transport, which requires a more complex treatment of heat transfer. This project implements a one-dimensional heat transfer solver based on the Boltzmann Transport Equation (BTE) and Monte Carlo methods to simulate phonon creation, annihilation, and transport. The solver is then extended to two dimensions with careful attention to the mesh Fourier number to improve accuracy.

## Key Concepts

- **Boltzmann Transport Equation (BTE)**: Describes the transport of phonons in materials.
- **Monte Carlo Simulations**: Used to model random processes such as phonon creation and annihilation events.
- **Phonon Transport**: Heat transfer in solids at the atomic scale, mediated by phonons.

## Features

- **1D Heat Transfer Solver**: Based on the BTE and Monte Carlo methods, simulating phonon transport.
- **2D Extension**: The solver is extended to 2D, incorporating geometric formulations for phonon propagation.
- **Mesh Fourier Number Consideration**: Ensures that the solver is accurate at the nanoscale by properly accounting for the mesh Fourier number.

## Requirements

- MATLAB
