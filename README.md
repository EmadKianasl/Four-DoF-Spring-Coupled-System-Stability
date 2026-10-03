# Stability Analysis of a Four-DOF Spring-Coupled Mechanical System

This repository contains the final project for the **Advanced Dynamics** course, focusing on the mathematical modeling, numerical simulation, and stability analysis of a nonlinear four-degree-of-freedom mechanical system consisting of four coupled rods with springs and damping.

The project develops the equations of motion using several formulations and investigates both the local and large-deviation stability characteristics of the system.

## Project Overview

The main objectives of this project are:

* Derivation of the nonlinear equations of motion
* Formulation using **Newton's method**
* Formulation using **Lagrange's equations**
* Application of the **Direct Hamilton Principle**
* Formulation using **canonical Hamiltonian equations**
* Numerical simulation of the dynamic response
* Linearization around equilibrium configurations
* Small-deviation stability analysis using eigenvalues
* Large-deviation stability analysis using **Lyapunov's method**
* Computation and visualization of Lyapunov exponent convergence
* Estimation and visualization of the system's region of attraction

## System Description

The studied system is a nonlinear mechanical system with four generalized coordinates:

$$
\boldsymbol{\theta}
=
[\theta_1,\theta_2,\theta_3,\theta_4]^T
$$

The four mechanical elements are coupled through elastic and damping effects. The model includes gravitational forces, spring coupling, viscous damping, and an external torque.

The governing equations are nonlinear due to the presence of the gravitational terms involving:

$$
\sin(\theta_i)
$$

The system is therefore suitable for studying both local linear stability and nonlinear stability characteristics.

## Stability Analysis

Two main stability approaches are investigated.

### 1. Small-Deviation Stability

The nonlinear equations are linearized around equilibrium configurations.

The resulting state-space model has the form:

$$
\dot{x}=Ax+Bu
$$

The eigenvalues of the system matrix are then calculated to determine the local stability characteristics of the equilibrium point.

The analysis considers equilibrium configurations around:

* \(\theta_i = 0\)
* \(\theta_i = \pi\)

For the selected parameters, the equilibrium around \(\theta_i=0\) is reported as locally asymptotically stable, while the configuration around \(\theta_i=\pi\) exhibits unstable small-deviation behavior.

## 2. Large-Deviation Stability

The nonlinear system is further investigated using a Lyapunov-based approach.

A quadratic Lyapunov function is constructed in the form:

$$
V(x)=x^TPx
$$

where \(P\) is a positive-definite matrix.

The time derivative of the Lyapunov function is investigated to determine the stability properties of the nonlinear system.

The project also examines the **region of attraction** associated with the selected Lyapunov function.

## Numerical Simulation

The dynamic response of the system is obtained numerically for a representative set of parameters and initial conditions.

The simulations include the time histories of:

* Angular displacement \(\theta_1\)
* Angular velocity \(\dot{\theta}_1\)
* Angular displacement \(\theta_2\)
* Angular velocity \(\dot{\theta}_2\)
* Angular displacement \(\theta_3\)
* Angular velocity \(\dot{\theta}_3\)
* Angular displacement \(\theta_4\)
* Angular velocity \(\dot{\theta}_4\)

The project also includes three-dimensional visualizations related to the Lyapunov analysis and region of attraction.

## Lyapunov Exponents

The convergence of the calculated Lyapunov exponent is investigated numerically.

The corresponding MATLAB implementation uses the system Jacobian and a quadratic Lyapunov formulation to investigate the stability characteristics of the nonlinear system.

## Project Structure

```text
four-dof-spring-coupled-system-stability/
│
├── README.md
│
├── MATLAB/
│   ├── equations_of_motion/
│   ├── simulation/
│   ├── linear_stability/
│   ├── lyapunov_analysis/
│   └── region_of_attraction/
│
├── figures/
│   ├── dynamic_response/
│   ├── eigenvalue_analysis/
│   └── lyapunov_analysis/
│
├── report/
│   └── Advanced_Dynamics_Project.pdf
│
└── LICENSE
```

## Tools

* **MATLAB**
* Symbolic and numerical computation
* Numerical integration of nonlinear equations
* Eigenvalue analysis
* State-space modeling
* Lyapunov stability analysis
* Data visualization

## Course Project

**Course:** Advanced Dynamics
**Department:** Mechanical Engineering
**Institution:** Isfahan University of Technology

**Student:** Emad Kian Asl
**Supervisor:** Dr. Mousavi, Associate Professor

## References

The theoretical development and numerical analysis are based on the course project report and the references cited therein.

## License

This project is intended primarily for academic and educational purposes.
