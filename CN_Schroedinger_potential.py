# -*- coding: utf-8 -*-
"""
Created on Sun Aug  7 19:53:25 2022

@author: KIERAN.WHITE

Crank-Nicolson method for 1D Schroedinger equation.
Working in natural units, with hbar = mass = 1.
Gaussian quantum particle with initial momentum kick p = m*v = hbar*k = k, 
where k is wave number. (a.k.a. Gaussian wave packet)
Lattice potential (but easily adjustable).
(row, col)
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Define function for norm.
def normf(dxc, uc, ic):
    return sum(dxc * np.square(np.abs(uc[ic, :])))

# Define function for expectation value of position.
def xexpf(dxc, xc, uc, ic):
    return sum(dxc * xc * np.square(np.abs(uc[ic, :])))

# Define function for expectation value of squared position.
def xexpsf(dxc, xc, uc, ic):
    return sum(dxc * np.square(xc) * np.square(np.abs(uc[ic, :])))

# Define function for standard deviation.
def sdaf(xexpc, xexpsc, ic):
    return np.sqrt(xexpsc[ic] - np.square(xexpc[ic]))

# Time t: t0 =< t =< tf. Have N steps at which to evaluate the CN scheme. The
# time interval is dt. decp: variable for plotting to certain number of decimal
# places.
t0 = 0
tf = 50
N = 500
dt = tf / N
t = np.linspace(t0, tf, num = N + 1, endpoint = True)
decp = str(dt)[::-1].find('.')

# Initialise array for filling with norm values at each time step.
norm = np.zeros(len(t))

# Initialise array for expectation value of position.
xexp = np.zeros(len(t))

# Initialise array for expectation value of squared position.
xexps = np.zeros(len(t))

# Initialise array for alternate standard deviation.
sda = np.zeros(len(t))

# Position x: -a =< x =< a. M is an even number. There are M + 1 total discrete
# positions, for the points to be symmetric and centred at x = 0.
a = 100
M = 1200
dx = (2 * a) / M
x = np.linspace(-a, a, num = M + 1, endpoint = True)

# The gaussian function u diffuses over time. sd sets the width of gaussian. u0
# is the initial gaussian at t0.
sd = 1
var = np.power(sd, 2)
mu = 0
k = 1
u0 = np.sqrt(1 / np.sqrt(np.pi * var)) * np.exp(-np.power(x - mu, 2) / (2 * \
     var)) * np.exp(1j * k * x)
u = np.zeros([len(t), len(x)], dtype = 'complex_')
u[0, :] = u0

# Normalise u.
u[0, :] = u[0, :] / np.sqrt(normf(dx, u, 0))

# Setup potential V. Want V(x < 0) = 0.
V0 = 1
l = 50 * sd
V = np.zeros(len(x))
for iv in range(len(x)):
    if np.sign(x[iv]) == 1:
        V[iv] = V0 * np.square(np.sin(2 * np.pi * x[iv] / l))
    elif np.sign(x[iv]) == -1:
        V[iv] = 0
    else:
        V[iv] = 0


# Set coefficients of CN scheme.
alpha = dt * -1j / (4 * np.power(dx, 2))
beta = dt * 1j * V / 2

# Tridiagonal matrices Al and AR. Al to be solved using Thomas algorithm.
Al = np.zeros([len(x), len(x)], dtype = 'complex_')

for i in range (0, M):
       Al[i + 1, i] = alpha
       Al[i, i] = 1 - (2 * alpha) + beta[i]
       Al[i, i + 1] = alpha

# Corner elements for BC's.
Al[M, M], Al[0, 0] = 1 - alpha + beta[i], 1 - alpha + beta[i]

Ar = np.zeros([len(x), len(x)], dtype = 'complex_')

for i in range (0, M):
       Ar[i + 1, i] = -alpha
       Ar[i, i] = 1 + (2 * alpha) - beta[i]
       Ar[i, i + 1] = -alpha

# Corner elements for BC's.
Ar[M, M], Ar[0, 0] = 1 + alpha - beta[i], 1 + alpha - beta[i]

# Thomas algorithm variables. Following similar naming as in Wikipedia page.
a = np.diag(Al, -1)
b = np.diag(Al)
c = np.diag(Al, 1)

NT = len(b)

cp = np.zeros(NT - 1, dtype = 'complex_')

for n in range(0, NT - 1):
    if n == 0:
        cp[n] = c[n] / b[n]
    else:
        cp[n] = c[n] / (b[n] - (a[n - 1] * cp[n - 1]))

d = np.zeros(NT, dtype = 'complex_')

dp = np.zeros(NT, dtype = 'complex_')

# Iterate over each time step to solve CN method. Maintain boundary
# conditions. Keep track of multiple values.
for i in range(0, N):
       # BC's.
       u[i, 0], u[i, M] = 0, 0
       
       # Find RHS.
       d = np.dot(Ar, u[i, :])
       
       for n in range(0, NT):
           if n == 0:
               dp[n] = d[n] / b[n]
           else:
               dp[n] = (d[n] - (a[n - 1] * dp[n - 1])) / (b[n] - (a[n - 1] * \
                                                                  cp[n - 1]))

       nc = NT - 1
       while nc > -1:
           if nc == NT - 1:
               u[i + 1, nc] = dp[nc]
               nc -= 1
           else:
               u[i + 1, nc] = dp[nc] - (cp[nc] * u[i + 1, nc + 1])
               nc -= 1

       norm[i] = normf(dx, u, i)
       xexp[i] = xexpf(dx, x, u, i)
       xexps[i] = xexpsf(dx, x, u, i)
       sda[i] = sdaf(xexp, xexps, i)


# Fill in final norm value.
norm[N] = normf(dx, u, N)

# Fill in final position expectation value.
xexp[N] = xexpf(dx, x, u, N)

# Fill in final squared position expectation value.
xexps[N] = xexpsf(dx, x, u, N)

# Fill in final standard deviation value.
sda[N] = sdaf(xexp, xexps, N)

# Animate the Gaussian evolution over time.
fig1, ax1 = plt.subplots(figsize = (8, 6))

def animatex(i):
       ax1.clear()
       ax1.set_ylim(-np.max(np.abs(u)), np.max(np.abs(u)))
       #ax1.plot(x, np.real(u[i, :]), color = 'r', label = 'real part')
       #ax1.plot(x, np.imag(u[i, :]), color = 'b', label = 'imagnary part')
       ax1.plot(x, np.abs(u[i, :]), color = 'm', label = 'absolute')
       ax1.plot(x, V, color='y')
       ax1.plot(0, 0, label = 'time t = {:.3f} seconds')#.format(round(t[i], \
       #decp)), color = 'w')
       ax1.set_xlabel(r'$x$')
       ax1.set_ylabel('amplitude')
       ax1.legend()

ani = FuncAnimation(fig1, animatex, repeat = True, interval=1e-24, frames = N \
                    + 1)

# Plot time evolution of distribution's standard deviation.
fig2, ax2 = plt.subplots(figsize = (8, 6))


ax2.plot(t, np.abs(sda), color = 'b', label = 'CN scheme 'r'$\sigma$''('r'$t$'\
         ')')
ax2.plot(t, sd * np.sqrt(1 + np.power(t / var, 2)) / np.sqrt(2), color = 'r', \
         label = 'analytical 'r'$\sigma$''('r'$t$'')')
ax2.set_title('Standard deviation 'r'$\sigma$'' as a function of time 'r'$t$')
ax2.set_xlabel(r'$t$'' (s)')
ax2.set_ylabel(r'$\sigma$')
ax2.legend()
