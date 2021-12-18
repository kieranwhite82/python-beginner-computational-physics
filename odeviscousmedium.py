"""
Program that solves an ODE using forward finite difference method.
A particular case: simulating a sphere falling through an infinte viscous
fluid. The forces taken into account are drag, weight and bouyancy.

Input
Initial velocity.
Number of intervals.
Final time.
Sphere radius.
Sphere mass.
Density of viscous fluid.

Outut
Sphere's velocity at each time interval.
Fluid characteristics: dynamic and kinematic viscosity.
"""
import numpy as np
import matplotlib.pyplot as plt

print("This program solves a first order ODE using forward finite difference m\
ethod. The sphere travels vertically. Upwards is positive, downwards negative.\
 There are several cases. If the sphere is more dense than the fluid, it shoul\
d reach a negative terminal velocity. If the sphere has the same density as th\
e fluid, it will reach zero velocity. If the sphere is less dense than the flu\
id, it will reach a postive terminal velocity.\n\nEnter all values in SI units\
.")

# Get user input for values.
v0 = eval(input("Enter the intial velocity: "))
N = eval(input("Enter the number of intervals: "))
tf = eval(input("Enter the final time: "))
sR = eval(input("Enter the radius of the sphere: "))
sM = eval(input("Enter the mass of the sphere: "))
fRho = eval(input("Enter the density of the fluid: "))

# Assign values to constants.
g = 9.81
# Sphere characteristics.
sV = (4 / 3) * np.pi * np.power(sR, 3)
sRho = sM / sV
# Drag coefficient.
Cd = 0.47
sA = np.pi * np.power(sR, 2)
b = fRho * Cd * sA / 2
# Time interval.
dt = tf / N
# Initialise arrays for storing time and velocity values.
v = np.zeros(N + 1)
v[0] = v0
t = np.zeros(N + 1)
t[0] = 0


# Forward finite difference function:
def ffd(i, dt):
       return v[i] + ((((sV * g / sM) * (fRho - sRho)) - (b * v[i] * abs(v[i])\
                      / sM)) * dt)


# Iterate to estimate velocity at each time interval.
i = 0
while i < N:
       v[i + 1] = ffd(i, dt)
       t[i + 1] = t[i] + dt
       i += 1


# Estimate fluid characteristics, dynamic and kinematic viscosity. Assuming
# terminal velocity has been reached by the end of the simulation.
fMu = (2 / 9) * (sRho - fRho) * g * np.power(sR, 2) / abs(v[N])
fNu = fMu / fRho

print("\nThe estimate of the dynamic viscosity: {0:.6g} PA.s\n\nThe estimate o\
f the kinematic viscosity: {1:.5g} J.s/kg".format(fMu, fNu))
# Plot derivative approximations.
plt.figure()
plt.xlabel("time (s)")
plt.ylabel("velocity (m/s)")
plt.plot(t, v, 'b', label='v(t)')
plt.legend()
