#!/usr/bin/env python3
from scipy.linalg import pascal
import mpmath as mp

# =============================================================================
# A-STABLE BDF NORDSIEK COEFFICIENTS
# =============================================================================
# Stability of the BDF method depends on the spectral radius of the matrix:
#
#    M(x) = P - P l (tau1 - x tau) / (r1 - x r),
#
#    where x = a h,
#          h                   : time step,
#          P                   : upper-triangular Pascal matrix,
#          l   = (l0,l1,...lq) : method's coefficients,
#          tau = (1,1,...,1)   : unit vector
#          tau1= (0,1,...,q)
#          r   = tau.l
#          r1  = tau1.l
#          q > 0               : method's order
#
# Here the spectrum of M is computed at x->inf for standard coefficients
# in an arbitrary precision.
# =============================================================================

mp.mp.dps = 90  # set precision (digits)

ll = {
    1 : mp.matrix ([0,  1]),
    2 : mp.matrix ([0,  1,  1]),
    3 : mp.matrix ([0,  2,  3,  1]),
    4 : mp.matrix ([0,  6, 11,  6,  1]),
    5 : mp.matrix ([0, 24, 50, 35, 10,  1]),
    6 : mp.matrix ([0,120,274,225, 85, 15, 1]),
}

for q in ll:

    P = pascal (q+1, kind='upper')

    P = mp.matrix (P)

    l = ll[q]

    e = mp.matrix([1 for i in range(q+1)])

    r = (e.T @ l)[0]

    u = P @ l

    M = P - u @ e.T / r

    res = mp.eig(M)[0]

    print (f"Order q={q:} : Re Im Abs")
    for i, a in enumerate(res):
        print(f"{i:}{a.real:12.3e}{a.imag:12.3e}{mp.fabs(a):12.3e}")
    print("")