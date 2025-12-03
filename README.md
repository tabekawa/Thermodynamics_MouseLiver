Thermodynamics_MouseLiver
====
GLEAM (Gibbs free energy of reaction Landscape Estimation from metabolome Assisted by variance-covariance Matrix) estimates Gibbs free energy change of reaction from metabolomic data with uncertainties and missing values. The code to perform analysis regarding some physiological characteristics affected by Gibbs free energy change of reaction, such as Flux Control Coefficient (FCC) and Enzyme Cost (EC) of the raction, are also provided.

# Requirement

The developmental version of package has been tested on the following systems:

Windows 11

MATLAB R2021b with Optimization Toolbox, Global Optimization Toolbox

Gurobi 10.0.3

# Set up

Download the package code to your local path.

# Contents

## GLEAM

This is the package code to estimate thermodynamically consistent metabolite concentrations, standard Gibbs free energy of formation, and Gibs free energy change of reaction by GLEAM.

## FCC

This is the package code to calculate FCC within glycolysis gluconeogenesis, and the TCA cycle.

## EC

This is the package code to calculate EC, Enzyme Cost Minimum (ECM), and Metabolite Equilibrium Gap (MEG) for glycolysis gluconeogenesis, and the TCA cycle.

## Mouse

The input data for the analysis performed in Abekawa et al. is in this directory.

# Demo

To perform the analysis in Abekawa et al., execute addpath(genpath('Mouse')) in the Thermodynamics_MouseLiver directory within MATLAB, and run demo_GLEAM_mouse.m, demo_FCC_mouse.m, and demo_EC_mouse.m in the GLEAM, FCC, and EC directory, respectively.

# Contact

Takumi Abekawa: [abekawa@g.ecc.u-tokyo.ac.jp](abekawa@g.ecc.u-tokyo.ac.jp)

Satoshi Ohno: [s-ohno.dsc@tmd.ac.jp](s-ohno.dsc@tmd.ac.jp)

Shinya Kuroda: [skuroda@bs.s.u-tokyo.ac.jp](skuroda@bs.s.u-tokyo.ac.jp)

# Reference
[1]User Name, 'Paper Titile' Conference Name pp.xx 20XX

[test](https://github.com/test)
