#!/usr/bin/env python3
import numpy as np
from ninteg import integrate

# =============================================================================
# SKINATH MODEL DESCRIPTION
# =============================================================================
#
# 1. Point Kinetics Equations (Neutron Power & Precursors):
#    dP/dt = (R(T) - 1) * [Σ(bet) / L] * P(t) + Σ(lam * C(t))
#    dC/dt = (bet / L) * P(t) - lam * C(t)
#
# 2. Thermal-Hydraulics Equation (Temperature):
#    C0 * dT/dt = P(t) - K0 * (1 - T0 / T(t))^(1/4) * (T(t) - T0)
#
# 3. Reactivity Feedback & Heat Transfer:
#    R(T) = R0 + B * (T - T0)
#    K0   = A0 * H0 * D0^(3/4)  (Convective Heat Transfer Coefficient)
#
# 4. Initial conditions: steady state at P0 and T0
# =============================================================================

# --- Physical Dimensions & Environmental Constants ---
H0 = 23.e-2          # Cylinder height (m)
D0 = 20.e-2          # Cylinder diameter (m)
T0 = 20.0            # Initial/Ambient temperature (°C)
A0 = 17.52           # Empirical physical constant for air

# --- Reactor Kinetics Parameters ---
P0 = 1.0E-2          # Initial Power (W)
L0 = 5.0E-5          # Prompt neutron lifetime (s)
R0 = 0.043           # Initial reactivity perturbation ($)
B0 = -3.06E-3        # Temperature reactivity coefficient ($/K)

# --- Derived Thermal Properties ---
V0 = np.pi * H0 * (D0**2) / 4     # Core volume (m^3)
C0 = 1.8e+6 * V0                  # Total heat capacity (Erg/°C)
K0 = A0 * H0 * (D0**0.75)         # Heat transfer coefficient

# --- Delayed Neutron Data (Precursor Groups) ---
lam = np.array([1.24e-2, 3.05e-2, 1.11e-1, 3.01e-1, 1.14e+0, 3.01e+0]) # Decay constants (s^-1)
bet = np.array([2.20e-4, 1.42e-3, 1.27e-3, 2.57e-3, 7.50e-4, 2.70e-4]) # Delayed neutron fractions

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
