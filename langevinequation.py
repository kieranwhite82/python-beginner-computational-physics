"""
Program that solves a second-order stochastic ordinary differential equation in
one dimension, for N particles with no potential acting on them.
"""
import random as rd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
"""
Initialise seed for random number generation, using default argument of time.
"""
rd.seed()
"""
Box-Mueller algorithm. Returns a randomly picked number from a Gaussian
distribution, of which the width is a function of variables dt and p. Slower in
comparison to calling ready-made function, but shows process of creating
pseudo-random samples with Gaussian distribution. Slow as every time the
function is called, a new distribution is created from which a number is picked
at random. This is the noise term; source of stochastic behaviour.
"""
def dW(N, sd):
       
       r1 = np.zeros(N)
       r2 = np.zeros(N)
       r1s = np.zeros(N)
       r2s = np.zeros(N)
       r = np.zeros(N)
       j = 0
       s = []
       x1 = []
       
       for i in range(0, N):
              r1[i] = rd.random()
              r2[i] = rd.random()
              r1s[i] = (2 * r1[i]) - 1.0
              r2s[i] = (2 * r2[i]) - 1.0
              r[i] = r1s[i]**2 + r2s[i]**2
       
              if r[i] > 0 and r[i] < 1:
                     S = np.sqrt(-2 * np.log(r[i]) / r[i])
                     s.append(S)
                     X1 = r1s[i] * s[j] * sd
                     x1.append(X1)
                     j += 1
              else:
                     pass

       return rd.choice(x1)

"""
Forward finite difference for velocity (Euler):
"""
def vl(i):
       return v[i] + (-l * v[i] * dt) + (g * dW(N, sd))

"""
Forward finite difference for position (Euler):
"""
def xl(i):
       return x[i] + (v[i] * dt)

"""
Variables for Langevin equation. Working in natural units: m, kB are taken to
be unity.
"""
l, T, m, kB = 10.0, 2.0, 1.0, 1.0
"""
Strength of noise term, g.
"""
g = np.sqrt(2 * l * kB * T / m)
"""
Final time to iterate to, tf. Initial time is treated as zero.
"""
tf = 10
"""
Number of time step intervals, N.
"""
N = 100
"""
Index for scaling term for determining width of Gaussian distribution of random
numbers for stochastic term, p.
"""
p = 1 / 2
"""
Time step size, dt.
"""
dt = tf / N
"""
Width of Gaussian term, sd.
"""
sd = np.power(dt, p)
"""
Initial x-position, x0.
"""
x0 = 0
"""
Initial velocity, v0.
"""
v0 = 0
"""
Number of times to simulate particle path, Nt.
"""
Nt = 250
"""
Number of intervals either side of mean value for both position and velocity to
bin data into, Nb. Three standard deviations, each split into four bins makes
twelve.
"""
Nb = 12
"""
Standard deviation of position, dx.
"""
dx = np.sqrt(2 * kB * T * tf / (l * m))
"""
Standard deviation of velocity, dv.
"""
dv = np.sqrt(kB * T / m)
"""
Initialise arrays, initial values.
"""
"""
Velocity, v.
"""
v = np.zeros(N + 1)
v[0] = v0
"""
Position, x.
"""
x = np.zeros(N + 1)
x[0] = x0
"""
Time, t.
"""
t = np.zeros(N + 1)
t[0] = 0

"""
Expectation value of position, Ex.
"""
Ex = np.zeros(N + 1)
Ex[0] = x[0]
"""
Expectation value of square of position, Ex2.
"""
Ex2 = np.zeros(N + 1)
Ex2[0] = np.power(x[0], 2)
"""
Variance of position, Vx.
"""
Vx = np.zeros(N + 1)
Vx[0] = Ex2[0] - np.square(Ex[0])

"""
Two-dimensional arrays for storing data.
"""
"""
Position storage.
"""
df = np.zeros([N + 1, Nt])
Pd = np.zeros([N + 1, 2 * Nb])
"""
Velocity storage.
"""
vf = np.zeros([N + 1, Nt])
Vd = np.zeros([N + 1, 2 * Nb])

"""
Solve Langevin equation numerically.
"""
for k in range(0, Nt):

       for i in range(0, N):
              v[i + 1] = vl(i)
              x[i + 1] = xl(i)
              t[i + 1] = (i + 1) * dt
              
              """
              Bin data given previous outlined variables at each time step.
              """
              for n1 in range(0, 2 * Nb):
                     
                     if x[i + 1] - x0 >= ((-Nb + n1) * dx / 4) and x[i + 1] - x0 < ((-Nb + n1 + 1) * dx / 4):
                            
                            Pd[i + 1, n1] += 1
                     else:
                            pass
              
              for n2 in range(0, 2 * Nb):
                     
                     if v[i + 1] - v0 >= ((-Nb + n2) * dv / 4) and v[i + 1] - v0 < ((-Nb + n2 + 1) * dv / 4):
                            
                            Vd[i + 1, n2] += 1
                     else:
                            pass
       """
       Store data for each trajectory.
       """
       df[:, k] = x
       vf[:, k] = v

"""
decp: variable for plotting to certain number of decimal places.
"""
decp = str(dt)[::-1].find('.')

"""
Px and lblx: lists for plotting position bars, given previous outline for binning data.
"""
Px = []
for i in range(0, (2 * Nb)):
       
       PX = dx * (-Nb + i) / 4
       Px.append(PX)


lblx = []
for i in range(0, int((Nb / 2) + 1)):
       if i == (Nb / 4):
              LBL = r"$\mu_x$"
              lblx.append(LBL)
       else:
              LBL = r"${0}\sigma_x$".format(-(Nb / 4) + i)
              lblx.append(LBL)

"""
Pv and lblv: lists for plotting velocity bars, given previous outline for binning data.
"""
Pv = []
for i in range(0, (2 * Nb)):
       
       PV = dv * (-Nb + i) / 4
       Pv.append(PV)


lblv = []
for i in range(0, int((Nb / 2) + 1)):
       if i == (Nb / 4):
              LBL = r"$\mu_v$"
              lblv.append(LBL)
       else:
              LBL = r"${0}\sigma_v$".format(-(Nb / 4) + i)
              lblv.append(LBL)

"""
Calculate velocity variance, Vv, given position expectation values.
"""
for i in range(0, N + 1):
       Ex[i] = np.mean(df[i, :])
       Ex2[i] = np.mean(np.power(df[i, :], 2))
       Vx[i] = Ex2[i] - np.power(Ex[i], 2)


Vv = l * Ex2[1:N + 1] / (2 * t[1:N + 1])


"""
Plotting figures. Static plots grouped, animated plots separate.
"""

fig1, axs = plt.subplots(2, 2)


for i in range(0, 2):
       
       for j in range(0, 2):
              
              axs[i, j].set_xlim(0, tf)


for k in range(0, Nt):
       axs[0, 0].plot(t, df[:, k])
axs[0, 0].set_title("Positions of {} trajectories".format(Nt))
axs[0, 0].set_ylabel(r"$x$")


for k in range(0, Nt):
       axs[0, 1].plot(t, vf[:, k])
axs[0, 1].set_title("Velocities of {} trajectories".format(Nt))
axs[0, 1].set_ylabel(r"$v$")


axs[1, 0].plot(t, Vx)
axs[1, 0].plot(t, 2 * kB * T * t / (l * m), label = "Theoretical "r"$\sigma_x$")
axs[1, 0].set_title("Position variance")
axs[1, 0].set_ylabel(r"$\sigma_x$")
axs[1, 0].set_xlabel(r"$t$"" (s)")
axs[1, 0].legend()


axs[1, 1].plot(t[1:N + 1], Vv)
axs[1, 1].plot(t, kB * T * np.ones(len(t)) / m, label = "Theoretical "r"$\sigma_v$")
axs[1, 1].set_title("Velocity variance")
axs[1, 1].set_ylabel(r"$\sigma_v$")
axs[1, 1].set_xlabel(r"$t$"" (s)")
axs[1, 1].legend()


fig1, ax1 = plt.subplots(figsize=(8, 6))

def animatex(i):
       ax1.clear()
       ax1.bar(Px, Pd[i, :], width=dx / 4, align='edge', color='b', \
               label = 't = {} seconds'.format(round(t[i], decp)))
       s_ticks = np.linspace(-3 * dx, 3 * dx, num = 7, endpoint = True)
       ax1.set_xticks(s_ticks)
       ax1.set_xticklabels(lblx)
       ax1.set_ylim(0, np.max(Pd))
       ax1.set_xlim(-3 * dx, 3 * dx)
       ax1.set_ylabel("Particle count")
       ax1.set_xlabel(r"$x$")
       ax1.legend()

anix = FuncAnimation(fig1, animatex, repeat = True, interval=200, frames=N + 1)


fig2, ax2 = plt.subplots(figsize=(8, 6))

def animatev(i):
       ax2.clear()
       ax2.bar(Pv, Vd[i, :], width=dv / 4, align='edge', color='b', \
               label = 't = {} seconds'.format(round(t[i], decp)))
       s_ticks = np.linspace(-3 * dv, 3 * dv, num = 7, endpoint = True)
       ax2.set_xticks(s_ticks)
       ax2.set_xticklabels(lblv)
       ax2.set_ylim(0, np.max(Vd))
       ax2.set_xlim(-3 * dv, 3 * dv)
       ax2.set_ylabel("Particle count")
       ax2.set_xlabel(r"$v$")
       ax2.legend()

aniv = FuncAnimation(fig2, animatev, repeat = True, interval=200, frames=N + 1)

plt.show()
