# Monte Carlo Simulations for Non-Fourier Heat Transfer Using Boltzmann Transport Equation


This repository contains the code and report for a research project focused on modeling heat transfer at the nanoscale, where Fourier's law of heat conduction no longer holds. The project extends Peraud's one-dimensional heat transfer solver based on the Boltzmann Transport Equation (BTE) and Monte Carlo methods to simulate phonon creation, annihilation, and transport. The solver is extended to two dimensions, with careful attention to the mesh Fourier number to improve accuracy.

## Overview

At the nanoscale, phonon transport deviates from the classical diffusive behavior described by Fourier's law. Instead, phonons exhibit ballistic or quasi-ballistic transport, which requires a more complex treatment of heat transfer. This project extends the one-dimensional heat transfer solver proposed by Peraud, based on the BTE and Monte Carlo methods, to simulate phonon transport. The extension to two dimensions incorporates geometric formulations for phonon propagation and ensures accuracy by considering the mesh Fourier number.

## Key Concepts

- **Boltzmann Transport Equation (BTE)**: Describes the transport of phonons in materials.
- **Monte Carlo Simulations**: Used to model random processes such as phonon creation and annihilation events.
- **Phonon Transport**: Heat transfer in solids at the atomic scale, mediated by phonons.

## Features

- **1D to 2D Extension**: Extends JP Peraud's 1D heat transfer solver to two dimensions for more accurate modeling of phonon transport.
- **Mesh Fourier Number Consideration**: Ensures that the solver is accurate at the nanoscale by properly accounting for the mesh Fourier number.

## Requirements

- MATLAB
