# NINTEG – Numerical Integration of Stiff ODEs

**NINTEG** is used to integrate stiff systems of Ordinary Differential Equations (ODEs) using Backward Differentiation Formula (BDF).
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
  the code handles nonlinear part iteratively.

## Usage
Solve dx / dt = a x on time interval 0 < t < 1 for initial value x(0) = 1 and parameter a = -1.5,

```python
from ninteg import integrate

# Define problem
time, x0, a = (0, 1), [1], -1.5

# Solve x - h * f (t,x) = b with tolerance e
def dynamics (h, t, b, x, e):
    x [:] = b [:] / (1 - a * h)

# Integrate
for t, x, info in integrate (time, x0, dynamics):
    print (t, x)
```

## Installation

```bash
git clone https://github.com/alcfb/ninteg.git
cd ninteg
pip install .
```

## Troubleshooting

- Error "too many trials at one time step": the in-step equation does not converge properly; improve the custom solver, and try smaller initial step size.
- Slow integration with a small step size: if the problem is nonlinear, convergence depends on the fixed-point iteration; try improving the quality of the custom solver in 'dynamics'.

## Reference

If you use this project, please cite:

> A. Cherezov, A. Vasiliev, H. Ferroukhi,
"Application of Backward Differential Formula and Anderson’s method for multigroup diffusion transient equation", Annals of Nuclear Energy, Volume 210, 2025, 110837, 
> https://doi.org/10.1016/j.anucene.2024.110837
