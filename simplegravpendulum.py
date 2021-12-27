"""
Program that solves a second order ODE governing the motion of a harmonic
pendulum. The only force acting is the restoration force; no friction.
Sphere on end of string. Initial value problem. Basic Euler and Euler-Cromer
method compared to show Euler method violating energy conservation principle.
May be potential to use more efficient data management (dataframes?) but for
the purpose of tracking the values everything is separate.

Input:
Initial angle from vertical.
Initial rate of change of angle with time.
Mass of sphere.
Length of string.
Final time to run simulation to.
Number of intervals.

Output:
Angle from vertical as function of time.
Energy of sphere as function of time.
"""
import numpy as np
import matplotlib.pyplot as plt


# Forward finite difference for rate of change of angle (Euler):
def fthd(i, dt):
       return thd[i] - (((g / l) * th[i]) * dt)


# Forward finite difference for rate of change of angle (Euler-Cromer):
def fthdc(i, dt):
       return thdc[i] - (((g / l) * thc[i]) * dt)


# Forward finite difference for angle (Euler):
def fth(i, dt):
       return th[i] + (thd[i] * dt)


# Forward finite difference for angle (Euler-Cromer):
def fthc(i, dt):
       return thc[i] + (thdc[i + 1] * dt)


# Kinetic energy (Euler):
def K(i):
       return (1 / 2) * m * np.power(l, 2) * np.power(thd[i], 2)


# Potential energy (Euler):
def P(i):
       return (1 / 2) * m * g * l * np.power(th[i], 2)


# Kinetic energy (Euler-Cromer):
def Kc(i):
       return (1 / 2) * m * np.power(l, 2) * np.power(thdc[i], 2)


# Potential energy (Euler-Cromer):
def Pc(i):
       return (1 / 2) * m * g * l * np.power(thc[i], 2)


# Kinetic energy (Euler-Cromer):
def Krk(i):
       return (1 / 2) * m * np.power(l, 2) * np.power(thdrk[i], 2)


# Potential energy (Euler-Cromer):
def Prk(i):
       return (1 / 2) * m * g * l * np.power(thrk[i], 2)


# Runge-Kutta method:
def RK(i, dt):
       a1 = - (g / l) * thrk[i]
       w1 = thdrk[i]
       a2 = - (g / l) * (thrk[i] + (w1 * dt / 2))
       w2 = thdrk[i] + (a1 * dt / 2)
       a3 = - (g / l) * (thrk[i] + (w2 * dt / 2))
       w3 = thdrk[i] + (a2 * dt / 2)
       a4 = - (g / l) * (thrk[i] + (w3 * dt))
       w4 = thdrk[i] + (a3 * dt)
       w = thdrk[i] + ((dt / 6) * (a1 + (2 * (a2 + a3)) + a4))
       th = thrk[i] + ((dt / 6) * (w1 + (2 * (w2 + w3)) + w4))
       return w, th


print("\nThis program solves the equation of motion for a simple gravity pendu\
lum. A sphere swings from a string. The string is massless, inextensible and d\
oesn't slack. The sphere is a point mass. The motion is in a two-dimensional a\
rc. There are no losses, so the sphere will oscillate with constant amplitude,\
 and total energy.\n\nEnter all values in SI units. Enter angles in degrees: t\
hey are converted into radians.")

# Get user input for values.
th0 = eval(input("Enter the intial angle: "))
thd0 = eval(input("Enter the initial rate of change of angle: "))
m = eval(input("Enter the mass of the sphere: "))
l = eval(input("Enter the length of the string: "))
tf = eval(input("Enter the final time to simulate to: "))
N = eval(input("Enter the number of intervals: "))

# Assign values to constants.
g = 9.81
# Time interval.
dt = tf / N
# Initialise arrays for storing time and velocity values.
i = 0
th = np.zeros(N + 1)
th[0] = th0
thc = np.zeros(N + 1)
thc[0] = th0
thrk = np.zeros(N + 1)
thrk[0] = th0
thd = np.zeros(N + 1)
thd[0] = thd0
thdc = np.zeros(N + 1)
thdc[0] = thd0
thdrk = np.zeros(N + 1)
thdrk[0] = thd0
t = np.zeros(N + 1)
t[0] = 0

KE = np.zeros(N + 1)
KE[0] = K(i)
PE = np.zeros(N + 1)
PE[0] = P(i)
ET = np.zeros(N + 1)
ET[0] = KE[0] + PE[0]

KEc = np.zeros(N + 1)
KEc[0] = Kc(i)
PEc = np.zeros(N + 1)
PEc[0] = Pc(i)
ETc = np.zeros(N + 1)
ETc[0] = KEc[0] + PEc[0]

KErk = np.zeros(N + 1)
KErk[0] = Krk(i)
PErk = np.zeros(N + 1)
PErk[0] = Prk(i)
ETrk = np.zeros(N + 1)
ETrk[0] = KErk[0] + PErk[0]


# Iterate to estimate angle at each time interval.
while i < N:
       thd[i + 1] = fthd(i, dt)
       
       thdc[i + 1] = fthdc(i, dt)
       
       th[i + 1] = fth(i, dt)
       
       thc[i + 1] = fthc(i, dt)
       
       thdrk[i + 1] , thrk[i + 1] = RK(i, dt)
       
       KE[i + 1] = K(i + 1)
       PE[i + 1] = P(i + 1)
       ET[i + 1] = KE[i + 1] + PE[i + 1]
       
       KEc[i + 1] = Kc(i + 1)
       PEc[i + 1] = Pc(i + 1)
       ETc[i + 1] = KEc[i + 1] + PEc[i + 1]
       
       KErk[i + 1] = Krk(i + 1)
       PErk[i + 1] = Prk(i + 1)
       ETrk[i + 1] = KErk[i + 1] + PErk[i + 1]
       
       t[i + 1] = t[i] + dt
       i += 1


# Plot derivative approximations.
fig, axs = plt.subplots(2, 3)
fig.add_subplot(frameon=False)
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Time (seconds)")

label = ["Theta (degrees)", "Energy (joules)"]


for l, ax in zip(label, axs):
    ax[0].set_ylabel(l)


axs[0, 0].plot(t, th, 'b')
axs[0, 0].set_title('Euler theta (t)')
axs[0, 1].plot(t, thc, 'b')
axs[0, 1].set_title('Euler-Cromer theta (t)')
axs[0, 2].plot(t, thrk, 'b')
axs[0, 2].set_title('Runge-Kutta theta (t)')
axs[1, 0].plot(t, ET, 'r')
axs[1, 0].set_title('Euler Energy (t)')
axs[1, 1].plot(t, ETc, 'r')
axs[1, 1].set_title('Euler-Cromer Energy (t)')
axs[1, 2].plot(t, ETrk, 'r')
axs[1, 2].set_title('Runge-Kutta Energy (t)')
