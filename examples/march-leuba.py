#!/usr/bin/env python3
import numpy as np
from ninteg import integrate

# =============================================================================
# MARCH-LEUBA MODEL OF A BOILING WATER REACTOR 
# PEACH BOTTOM TEST CASE 3PT3
# =============================================================================
#
# 1. Point Kinetics Equations (Neutron Power & Precursors):
#    dP/dt = (rho - bet) / l0 * P + lam * C
#    dC/dt = bet / l0 * P - lam * C
#
# 2. Thermal-Hydraulics Equation (Temperature):
#    dT/dt = a1 * (P - 1) - a2 * T
#    dR/dt = Z
#    dZ/dt =-6 / tau Z - 12/tau^2 R + C H^2 / tau (dT/dt + 6/tau T)
#
# 3. Reactivity Feedback:
#    rho = R + D * T
#
# 4. Initial conditions: steady state at P0 = 1 and T0 = 0
# =============================================================================

class Point_Kinetics:

    l0  = 4.0E-5    # Prompt neutron lifetime (s)
    lam = [0.08]    # Decay constants (1/s)
    bet = [0.0056]  # Delayed neutron fractions (-)
    D0  =-2.61e-5   # Doppler Effect (1/K)
    R0  = 5.6e-5    # Reactivity Perturbation (-)

    def __init__ (self):
        self.bet = np.array (self.bet)
        self.lam = np.array (self.lam)
        rank = len(self.lam) + 1
        mat = np.zeros((rank, rank))
        mat[0,  0] =-sum(self.bet) / self.l0
        mat[0, 1:] = self.lam
        mat[1:, 0] = self.bet / self.l0
        np.fill_diagonal (mat[1:, 1:], -self.lam)
        self.mat_0 = mat
        self.mat_h = np.zeros((rank, rank))
        self.mat_e = np.eye(rank)
        self.src_h = np.zeros(rank)

    def solve (self, source, h, temp=0, void=0):
        self.mat_h [:,:] = self.mat_e - h * self.mat_0
        self.mat_h [0,0]-= h / self.l0 * (void + self.D0 * temp + self.R0)
        self.src_h [:] = source
        return np.linalg.solve (self.mat_h, self.src_h)

    def steady (self):
        return np.concatenate (([1], self.bet / self.lam / self.l0))

class Thermal_Hydraulics:

    a1 = 19.08      # K/s
    a2 =  0.19      # fuel temp. response delay (1/s)
    t0 =  1.63      # bubble transit time (s)
    C0 = -3.65e-4   # 1/K
    k0 =  1.7       # feedback gain (K/s)
    H0 =  4.        # height (m)

    def __init__ (self, rank=3):
        mat = np.zeros ((rank, rank))
        mat [0,0] =-self.a2
        mat [1,2] = 1
        mat [2,:] = [self.k0 * (6 / self.t0 - self.a2) * self.C0 * self.H0**2 / self.t0, -6 / self.t0, -12 / self.t0**2]
        self.mat_0 = mat
        self.mat_h = np.zeros ((rank, rank))
        self.mat_e = np.eye(rank)
        self.src_h = np.zeros(rank)
    def solve (self, source, h, power=0):
        self.mat_h [:,:] = self.mat_e - h * self.mat_0
        self.src_h [:] = source
        self.src_h [0] += h * self.a1 * (power - 1)
        self.src_h [2] += h * self.a1 * (power - 1) * self.C0 * self.H0**2 / self.t0
        return np.linalg.solve (self.mat_h, self.src_h)

    def steady (self):
        return np.zeros(3)


# Initialize Elementary Models
PK = Point_Kinetics ()
TH = Thermal_Hydraulics ()

def dynamics(h, t, b, x, e, params):
    x[:-3] = PK.solve (b[:-3], h, temp=x[-3], void=x[-2])
    x[-3:] = TH.solve (b[-3:], h, power=x[0])

init_state = np.concatenate ((PK.steady(), TH.steady()))

#solution = integrate((0, 60), init_state, dynamics, qmax=2, h0=1.E-2, control=False)

solution = integrate((0, 60), init_state, dynamics, rtol=1.E-3)

pcm = 1.E+5

for t, x, info in solution:
    print (f"t={t:9.3f}s : h={info.last_step_size:9.3e} s   q={info.last_order:1.0f}   P={x[0]:15.9e} W  RFUE={PK.D0 * x[-3] * pcm:5.1f} [pcm]  RVOI={x[-2]*pcm:5.1f} [pcm]  RTOT={(PK.R0 + PK.D0 * x[-3] + x[-2]) * pcm:5.1f} [pcm]")

# Reference solution
t_ref = 60
P_ref = 9.892487713e-01
T_ref = 1.105045744e-02

# Output Results
print(f"""   
          Time: {t:9.3f} (s)
         Power: {x[0] :19.9e} (-)    dP = {100*(1-x[0] /P_ref):6.1e} %
   Temperature: {x[-3]:19.9e} (K)    dT = {100*(1-x[-3]/T_ref):6.1e} %
   Total steps: {info.successful_steps:5.0f}
Rejected steps: {info.rejected_steps:5.0f}
Function calls: {info.function_calls:5.0f}""")
