#!/usr/bin/env python3
from ninteg import IVP

def dynamics (h, t, y, x, e):
    a = [1.0, -1.5]
    x[0] = (y[0] + h * x[1]**2) / (1 - a[0] * h)
    x[1] = (y[1] - h * x[0]**2) / (1 - a[1] * h)

x0 = [1, 1]

ivp = IVP ()

ivp.setup (0, x0, dynamics)

for t, x in ivp.solve (t=1.2, h=1.E-6):
    pass

print (f"Solution: T = {t:12.6e}, X = ({x[0]:12.6e}, {x[1]:12.6e})")
ivp.info()