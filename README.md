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

The resulting model is represented as an **8-state nonlinear dynamical system**, consisting of four angular positions and four angular velocities.

## Mathematical Formulations

The equations of motion are investigated using several formulations:

### Newton's Method

The equations of motion are derived from force and moment balances applied to the mechanical elements.

### Lagrange's Equations

The system is formulated using the Lagrangian:

$$
L=T-V
$$

and the corresponding Lagrange equations.

### Direct Hamilton Principle

The equations of motion are also developed using the variational formulation based on Hamilton's principle.

### Canonical Hamiltonian Formulation

The system is represented using canonical coordinates and momenta to obtain the corresponding Hamiltonian equations.

## Dynamic Response

The nonlinear equations are numerically integrated to investigate the dynamic response of the system for selected initial conditions and parameters.

The simulations include the time histories of:

* Angular displacement \(\theta_1\)
* Angular velocity \(\dot{\theta}_1\)
* Angular displacement \(\theta_2\)
* Angular velocity \(\dot{\theta}_2\)
* Angular displacement \(\theta_3\)
* Angular velocity \(\dot{\theta}_3\)
* Angular displacement \(\theta_4\)
* Angular velocity \(\dot{\theta}_4\)

The main nonlinear equations of motion used in the numerical simulations are implemented in `eom18.m`.

## Stability Analysis

The project investigates stability from two complementary perspectives.

### Small-Deviation Stability

The nonlinear equations are linearized around equilibrium configurations.

The resulting state-space model has the form:

$$
\dot{x}=Ax+Bu
$$

The eigenvalues of the system matrix are calculated to characterize the local stability of the equilibrium configurations.

The analysis considers equilibrium configurations around:

* \(\theta_i=0\)
* \(\theta_i=\pi\)

For the selected parameters, the equilibrium around \(\theta_i=0\) is reported as locally asymptotically stable, while the configuration around \(\theta_i=\pi\) exhibits unstable small-deviation behavior.

### Large-Deviation Stability

The nonlinear system is further investigated using a Lyapunov-based approach.

A quadratic Lyapunov function is constructed in the form:

$$
V(x)=x^TPx
$$

where \(P\) is a positive-definite matrix.

The time derivative of the Lyapunov function is analyzed to investigate nonlinear stability and the corresponding region of attraction.

## Lyapunov Exponents

The convergence of the Lyapunov exponent is investigated numerically.

The repository includes the **LET (Lyapunov Exponents Toolbox)** implementation used in the project. The analysis is based on the system Jacobian and Lyapunov-related matrix formulations to investigate the stability characteristics of the nonlinear system.

## Repository Structure

```text
Four-DoF-Spring-Coupled-System-Stability/
│
├── LET/
│   └── Lyapunov Exponents Toolbox
│
├── canonical hamilton/
│   └── Canonical Hamiltonian formulation
│
├── lagrange(initial1)/
│   └── Lagrangian formulation – initial condition set 1
│
├── lagrange(initial2)/
│   └── Lagrangian formulation – initial condition set 2
│
├── largperturbation_lyapunov/
│   └── Large-deviation and Lyapunov stability analysis
│
├── smallperturbation/
│   └── Small-deviation linear stability analysis
│
├── eom18.m
│   └── Nonlinear equations of motion in first-order state-space form
│
├── Four_Bar_Parallel_System_Stability_Analysis.pdf
│   └── Complete project report
│
└── README.md
```

## Main Components

| Directory / File                                          | Description                                                                             |
| --------------------------------------------------------- | --------------------------------------------------------------------------------------- |
| `eom18.m`                                                 | Nonlinear equations of motion represented as an 8-state first-order system              |
| `lagrange(initial1)`                                      | Numerical analysis based on the Lagrangian formulation for the first initial condition  |
| `lagrange(initial2)`                                      | Numerical analysis based on the Lagrangian formulation for the second initial condition |
| `canonical hamilton`                                      | Canonical Hamiltonian formulation and analysis                                          |
| `smallperturbation`                                       | Linearization and small-deviation stability analysis                                    |
| `largperturbation_lyapunov`                               | Large-deviation stability and Lyapunov analysis                                         |
| `LET`                                                     | Lyapunov exponent calculation and convergence analysis                                  |
| `Four_Bar_Parallel_System_Stability_Analysis.pdf` | Complete project report                                                         |

## Computational Tools

The project is implemented primarily in **MATLAB** and uses numerical computation for:

* Nonlinear ordinary differential equations
* Numerical integration
* State-space modeling
* Eigenvalue analysis
* Linearization
* Lyapunov stability analysis
* Lyapunov exponent calculation
* Dynamic-response visualization

## How to Use

### Requirements

* MATLAB
* Appropriate MATLAB toolboxes required by the individual scripts
* LET toolbox for the Lyapunov exponent analysis

### Basic Workflow

1. Clone the repository:

```bash
git clone https://github.com/EmadKianasl/Four-DoF-Spring-Coupled-System-Stability.git
```

2. Open the repository in MATLAB.

3. Start with `eom18.m` to inspect the nonlinear state-space representation of the system.

4. Explore the formulation-specific directories:

```text
lagrange(initial1)
lagrange(initial2)
canonical hamilton
```

5. For local stability analysis, use:

```text
smallperturbation
```

6. For nonlinear/Lyapunov stability analysis, use:

```text
largperturbation_lyapunov
```

7. For Lyapunov exponent calculations, use:

```text
LET
```

## Project Results

The project investigates:

* Nonlinear dynamic response of the four-DOF system
* Local stability around equilibrium configurations
* Eigenvalue-based stability characteristics
* Lyapunov-based nonlinear stability
* Convergence of Lyapunov exponents
* Region of attraction associated with the Lyapunov function

The complete derivations, mathematical development, numerical results, and figures are provided in the project report:

**`Four_Bar_Parallel_System_Stability_Analysis.pdf`**

## Course Project

**Course:** Advanced Dynamics
**Department:** Mechanical Engineering
**Institution:** Isfahan University of Technology

**Student:** Emad Kian Asl
**Supervisor:** Dr. Mousavi, Associate Professor

## References

The theoretical development and numerical analysis are based on the course project report and the references cited therein.

## License

This repository is intended primarily for academic and educational purposes.
