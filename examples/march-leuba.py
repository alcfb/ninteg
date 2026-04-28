#!/usr/bin/env python3
import numpy as np
from ninteg import integrate

# =============================================================================
# MARCH-LEUBA MODEL OF A BOILING WATER REACTOR 
# PEACH BOTTOM TEST CASE 3PT3
# =============================================================================
#
# 1. Point Kinetics Equations (Neutron Power & Precursors):
#    dP/dt = (R(T) - 1) * bet / L * P(t) + lam * C(t)
#    dC/dt = (bet / L) * P(t) - lam * C(t)
#
# 2. Thermal-Hydraulics Equation (Temperature):
#    dTf/dt = a1 * P(t) - a2 * (Tf(t) - T0)
#    dRa/dt + 6 / tau dRa/dt + 12/tau^2 Ra = C H^2 / tau (dQ/dt + 6/tau Q)
#
# 3. Reactivity Feedback:
#    R = Ra + D*Tf
#
# 4. Initial conditions: steady state at P0 and T0
# =============================================================================

#H = 4.D0
#C = -3.65D-4
#D = -2.61D-5
#tau = 1.63D0
#a3 = 6.D0 / tau
#a4 = 12.D0 / tau**2
#a5 = a1 * C * H**2 / tau
#k0 = C*H**2/tau * (a3-a2)
#k0 = -0.01798854 * 0.2

# --- Physical Dimensions & Environmental Constants ---
H = 4.          # Height (m)

# --- Reactor Kinetics Parameters ---
P0 = 1.0E-2     # Initial Power (W)
L0 = 4.0E-5     # Prompt neutron lifetime (s)
R0 = 0.01      # Initial reactivity perturbation ($)
D0 =-2.61E-5    # Doppler coefficient (1/K)

# --- Delayed Neutron Precursors ---
lam = [0.08]    # Decay constants (s^-1)
bet = [0.0056]  # Delayed neutron fractions (-)

# --- Thermal Hydraulic Properties ---
a1 = 19.08      # K/s
a2 =  0.19      # 1/s
dt =  1.63      # bubble transit time (s)
c0 = -3.65e-4   # 1/K
k0 = -0.0037    # K/s

# --- Computational Pre-calculations ---
sum_bet_L0 = np.sum(bet) / L0
rank = len(lam) + 1
E0 = np.eye(rank)

# Initialize Point Kinetics Matrix (PK)
PK = np.zeros((rank, rank))
PK[0, 0] = (R0 - 1) * sum_bet_L0
PK[0, 1:] = lam
PK[1:, 0] = bet / L0
np.fill_diagonal(PK[1:, 1:], -lam)

def dynamics(h, t, b, x, e):
    # Point Kinetics with Temperature Feedback
    Ah = E0 - h * PK
    Ah[0, 0] -= h * B0 * (x[-1] - T0) * sum_bet_L0
    x[:-1] = np.linalg.solve(Ah, b[:-1])

    # Thermal Hydraulics:
    C1 = K0 * h * (1 - T0 / x[-1])**0.25
    x[-1] = T0 + (C0 * (b[-1] - T0) + h * x[0]) / (C0 + C1)

t_ref = [10, 20, 10] # ?
p_ref = [-0.25115159e-01, 0.3530984878e+00, 0.72166432e-03]

# Initial state vector: [Power, Precursors, Temperature]
x0 = np.concatenate(([P0], P0 * bet / (lam * L0), [T0]))
solution = integrate((0, 250*60), x0, dynamics)

for t, x, info in solution: pass

# Output Results
print(f"""   
          Time: {t:9.3f} s
     Reference: 2.118680532E-03 Watt
         Power: {x[0]:15.9e} Watt
   Temperature: {x[-1]:9.3e} °C
   Total steps: {info.successful_steps:5.0f}
Rejected steps: {info.rejected_steps:5.0f}
Function calls: {info.function_calls:5.0f}""")
