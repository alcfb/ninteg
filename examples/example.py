#!/usr/bin/env python3
from ninteg import integrate

# Glycolytic oscillator
#
# dx/dt = - x + a y + x^2 y
# dy/dt =   b - a y - x^2 y
#
# Stable limit cycle: a = 0.08, b = 0.6
# Stable fixed point: a = 0.08, b = 1.0

def dynamics (h, t, b, x, e):
    A, B = 0.08, 1.0
    x[0] = (b[0] + A * h * x[1]) / (1 + h - h * x[0] * x[1])
    x[1] = (b[1] + B * h - h * x[0]**2 * x[1]) / (1 + A * h)

x0 = [1, 1]
time = (0, 100)

solution = integrate (time, x0, dynamics)

for t, x, info in solution: pass

print(f"""   
        Time T: {t:9.3f} s
    Solution X: ({x[0]:9.3e}, {x[1]:9.3e})
   Total steps: {info.successful_steps:5.0f}
Rejected steps: {info.rejected_steps:5.0f}
Function calls: {info.function_calls:5.0f}""")
