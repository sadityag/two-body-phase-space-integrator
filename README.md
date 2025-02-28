# two body phase space integrator
 FeynCalc Two-Body Integrator Module – A Mathematica module for performing two-body phase space integrations using FeynCalc. It provides pre-defined functions and a flexible module to integrate matrix elements squared, verified for FeynCalc 10 and cross‐checked against “From Spinors to Supersymmetry” Appendix D.
 # FeynCalc Two-Body Integrator Module

This repository contains a Mathematica module designed for performing two-body phase space integrations in FeynCalc. It provides a set of pre-defined functions for computing phase space factors as well as a flexible module that integrates a given matrix element squared over the two-body phase space. The formulae have been cross-checked against Appendix D of *From Spinors to Supersymmetry* (coauthored by Howard Haber).

> **Note:** This module has been verified to work with FeynCalc 10. Feel free to use it without any accreditation.

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Example](#example)
- [Contributing](#contributing)
- [License](#license)
- [References](#references)
- [Contact](#contact)

---

## Overview

The module facilitates the evaluation of phase space integrals for two-body decays or scatterings. It assumes that the matrix element squared is at most quadratic in the momenta of the outgoing particles and is defined in 3+1 dimensions (with a possible update to d dimensions in the future).

The main routine, `TwoBodyIntegratorModule`, takes as input:
- The matrix element squared.
- Masses for the two outgoing states (default is 0, but symbolic masses are allowed).
- Labels for the particle momenta (default: `p1` and `p2`).
- A four-vector (default: `q`) representing momentum conservation.

The integration is performed over the Lorentz-invariant two-body phase space with the delta function enforcing momentum conservation.

---

## Features

- **Pre-defined Functions:**  
  - `Kallenλ[aa, bb, cc]` for computing the Källén function.
  - Phase space integration functions like `twobody$R2`, `twobody$T1qμ`, `twobody$T11metric`, etc.

- **TwoBodyIntegratorModule:**  
  A versatile module that combines the above functions with derivatives (via FeynCalc’s `FourDivergence`) to perform the full phase space integration.

- **Example Application:**  
  An example demonstrating the evaluation of a sample matrix element, spin summation, and further simplification is provided in the code.

---

## Installation

### Prerequisites

- [Mathematica](https://www.wolfram.com/mathematica/)
- [FeynCalc](http://www.feyncalc.org/) (version 10 is recommended)

### Setup

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/yourusername/FeynCalc-TwoBodyIntegratorModule.git
   cd FeynCalc-TwoBodyIntegratorModule

