# NINTEG – Numerical Integration of Stiff ODEs

A stiff system of ODEs is integrated by Backward Differentiation Formula (BDF) method using automatic step size and order selection. The key feature of the code is that it allows users to provide their own in-step solver instead of a Jacobian matrix, giving more control over the integration process. In addition, only a linear part of this equation can be resolved, leaving the code to handle a nonlinear part iteratively using the Anderson method. The code does not require any external dependencies.

## Description

The code implements the BDF method for the system of ODEs, `dx/dt=f(t,x)`, where `x` is the unknown vector, `t` is the time, and `h` is the step-size. BDF method is an `A`-stable implicit multistep algorithm. In this implementation, solution is represented by the Nordsieck vector, `z = (x,hx',h^2x''/2!,...,h^q x^(q)/q!)`, allowing efficiently change step size `h` and order `q` of the integration scheme to control the local truncation error. The algorithm generally follows the LSODE method with that difference that the solution of in-step equation, `x - h f(t,x) = b`, is user-provided. The Anderson method is built-in to solve the in-step non-linear equation.

## Usage
Solve `dx/dt = a x` on time interval `0 < t < 1` for initial value `x(0) = 1` and parameter `a = -1.5`,

```python
from ninteg import integrate

# Define problem
time, x0, a = (0.0, 1.0), [1], -1.5

def dynamics (h, t, b, x, e):
    """
    User-defined implicit step function
    x - h * f(t, x) = b

    Parameters
    ----------
    h : float
        Current step size.
    t : float
        Current time.
    b : ndarray
        Right-hand side vector (given by solver).
    x : ndarray
        Solution vector (must be updated in-place).
        On entry, contains the current state.
        On exit, must contain the updated solution.
    e : float
        Requested tolerance
    """
    x [:] = b [:] / (1 - a * h)

# Integrate
for t, x, info in integrate (time, x0, dynamics):
    """
    integrate parameters
    --------------------
    time : (t0, t1)
        Time interval.
    x0 : array-like, shape `(n,)`
        Initial vector.
    dynamics : callable
        Function implementing implicit step.
    rtol : float, optional
        Relative tolerance, default 1.E-6.
    atol : float, optional
        Relative tolerance, default 1.E-10.
    h0 : float, optional
        Initial step size, default 1.E-6.
    hmin, hmax : float, optional
        Step size bounds, default from 1.E-10 to 1.E+10.
    qmax : int, optional
        Maximum method order, default 5.
    --------------------
    Results:
    t : float
        current time.
    x : array-like, shape `(n,)`
        State vector.
    info : Info
        Solver diagnostics container.
    """
    print (t, x)
```

## Installation

```bash
git clone https://github.com/alcfb/ninteg.git
cd ninteg
pip install .
```

## Troubleshooting

**If in-step iteration does not converge**
- May result in slow integration with very small time steps `h` and orders `q = 1`
- Furthermore, integration may fail with an error "too many trials at one time step"

**Fix:**
- Increase relative tolerance `rtol`
- Start with a smaller initial step size `h0`
- Improve convergence of in-step function `dynamics` (most often)
- Ensure that `f(t,x)` is a continuous function
- Apply preconditioning procedures

## Reference

If you use this project, please cite:

> A. Cherezov, A. Vasiliev, H. Ferroukhi,
"Application of Backward Differential Formula and Anderson’s method for multigroup diffusion transient equation", Annals of Nuclear Energy, Volume 210, 2025, 110837, 
> https://doi.org/10.1016/j.anucene.2024.110837
