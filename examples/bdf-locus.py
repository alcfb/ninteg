#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import pascal

ll = {
    1 : np.array ([0,  1]),
    2 : np.array ([0,  1,  1]),
    3 : np.array ([0,  2,  3,  1]),
    4 : np.array ([0,  6, 11,  6,  1]),
    5 : np.array ([0, 24, 50, 35, 10,  1]),
    6 : np.array ([0,120,274,225, 85, 15, 1]),
}

class StabRegion:
    def __init__ (self, l):
        """Construct M = P (I - l e^T)."""
        n = l.shape[0]
        self.P = pascal(n, kind='upper')
        self.I = np.eye(n)
        self.e = np.ones(n)
        self.l = np.array(l)
        self.u = self.P @ l
        self.u = self.u.reshape(n,1)
        self.de = np.arange(n)
        self.r = sum(l)
        self.dr = self.de @ self.l
        self.n = n

    def solve (self, x):
        s = (self.de - x * self.e) / (self.dr - x * self.r)
        s = s.reshape (self.n, 1)
        M = self.P - self.u @ s.T
        res = np.linalg.eigvals(M)
        return res


def plot_eigenvalues (eigvals):

    masked_eigvals = np.ma.masked_where (eigvals > 1, eigvals)

    fig, ax = plt.subplots (figsize=(5,5), dpi=160)

    ax.set_xlabel("Re")
    ax.set_ylabel("Im")

    ax.imshow (masked_eigvals, cmap='terrain')
    #ax.scatter(points.real, points.imag, s=1, c=eigvals, alpha=1, cmap='jet')

    #ax.grid(True)
    plt.tight_layout()
    plt.show()

# ---- Example usage ----
q = 6

l = np.random.rand(q+1)

l = ll [q]

stab = StabRegion (l)

k = 400
w = 120
X = np.linspace(-w, w, k)
Y = np.linspace(-w, w, k)

points = []
eigvals = []

for x in X:
    for y in Y:
        z = x + 1j * y
        eigvals = np.append (eigvals, max(np.abs(stab.solve(z))))
        points.append (z)

eigvals = np.array (eigvals).reshape (k, k).T

plot_eigenvalues (eigvals)


