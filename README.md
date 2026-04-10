# NINTEG – Numerical Integration of Stiff ODEs

## Description
A stiff system of ODEs, `dx/dt=f(t,x)`, is integrated step-by-step by Backward Differentiation Formula (BDF) in Nordsieck formulation, with automatic step size and order selection. According to this method, a non-linear system of algebraic equations, `x - h f(t,x) = b`, must be solved for every time step, where `x` is the unknown vector, `t` is the time, `h` is the step-size, and `b` is the source vector. The key feature of the code is that it allows users to provide their own in-step solver instead of a Jacobian matrix, giving more control over the integration process. Alternatively, only a linear part of this equation can be resolved, leaving the code to handle a nonlinear part iteratively using the Anderson method.

## Usage
Solve `dx/dt = a x` on time interval `0 < t < 1` for initial value `x(0) = 1` and parameter `a = -1.5`,

```python
from ninteg import integrate

# Define problem
time, x0, a = (0, 1), [1], -1.5

# Solve x - h * f (t,x) = b
# given step size h, time t and tolerance e
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
