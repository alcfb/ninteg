# NINTEG – Numerical Integrator for Stiff ODEs

**NINTEG** is to integrate stiff systems of ordinary differential equations (ODEs) using Backward Differentiation Formula (BDF).
The key feature is that it allows users to provide their own solver instead of a Jacobian matrix, giving more control over the integration process.

## Features

- **Stiff ODEs**  
  Stable and efficient integration for stiff systems, dx/dt=f(t,x), using BDF method in Nordsieck formulation.
  Automatically adjusts both step size and method order according to the LSODE algorithm.

- **Custom In-Step Solver**  
  Users supply their own solver for solution of in-step system: x - h f(t,x) = b,
  where x is the solution vector, t is the time, h is the step-size and b is the source vector.

- **Built-in Anderson Method**  
  The non-linear in-step equation is iterated using the Anderson method.
  Thus, only a linear part of the equation can be solved by the user, while
  the code handles nonlinear part automatically.

## Usage
Solve dx/dt = a x, a = -1.5, x(0) = 1, 0 < t < 1

```
from ninteg import integrate

# Define problem
time, x0, a = (0, 1), [1], -1.5

# Solve x - h * f (t,x) = b with tolerance e
def dynamics (h, t, b, x, e):
    x [:] = b [:] / (1 - a * h)

solution = integrate (time, x0, dynamics)

# Print solution
for t, x, info in solution:
    print (t, x)
```

## Installation

```bash
git clone https://github.com/alcfb/ninteg.git
cd ninteg
pip install .
```

## Tech Stack

- **Core:** Fortran  
- **Interface:** Python (wrapper)


## Reference

If you use this project, please cite:

> A. Cherezov, A. Vasiliev, H. Ferroukhi,
"Application of Backward Differential Formula and Anderson’s method for multigroup diffusion transient equation", Annals of Nuclear Energy, Volume 210, 2025, 110837, 
> https://doi.org/10.1016/j.anucene.2024.110837
