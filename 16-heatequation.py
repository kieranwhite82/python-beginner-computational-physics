"""
Program that solves the diffusion PDE through the FTCS method with the initial
condition of a gaussian. There is some uniform background concentration.
One-dimensional space. The diffusion coefficient is constant, so the PDE is 
identical to the heat equation. Thus concentration here is temperature.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Define function for calculating standard deviation of temperature 
# distribution.
def sdf(dxc, xc, uc, ic, ubackc):
       return np.sqrt(sum(dxc * np.power(xc, 2) * (u[ic, :] - (ubackc * \
       np.ones(len(xc))))) / sum(dxc * (u[ic, :] - (ubackc * \
       np.ones(len(xc))))))

# Time t: t0 =< t =< tf. Have N steps at which to evaluate the FTCS scheme. The
# time interval is dt. decp: variable for plotting to certain number of decimal
# places.
t0 = 0
tf = 5
N = 200
dt = tf / N
t = np.linspace(t0, tf, num = N + 1, endpoint = True)
decp = str(dt)[::-1].find('.')

# Position x: -a =< x =< a. M is an even number. There are M + 1 total discrete
# positions, for the points to be symmetric and centred at x = 0.
a = 20
M = 100
dx = (2 * a) / M
x = np.linspace(-a, a, num = M + 1, endpoint = True)

# The gaussian function u(x, t) diffuses over time. uback is the uniform
# background temperature. sd sets the width of gaussian. u0 is the initial
# gaussian at t0. upeak is the peak size of the inital temperature
# distribution.
uback = 1
upeak = 1.5
sd = 1
var = np.power(sd, 2)
u0 = upeak * np.exp(-np.power(x, 2) / (2 * var)) / np.sqrt(2 * np.pi * var)
u = np.zeros([len(t), len(x)])
u[0, :] = u0 + uback

# Set the diffusion constant, D.
D = 3.2

# From von Neumann stability analysis; the variable S is used to check if the
# solution is numerically stable.
S = D * dt / np.power(dx, 2)

print("\nThe stability variable is {}".format(S))

if S == 0.5:
       print("\nThe stability variable is equal to 1/2, shape should be \
representative of true solution and solution should be stable.")
elif S < 1:
       print("\nThe stability variable is less than 1/2, shape will not be \
representative of true solution and solution will be stable.")
else:
       print("\nThe stability variable is greater than 1/2, shape will not be \
representative of true solution and solution will not be stable.")

# Tridiagonal matrix A to solve set of linear equations making up FTCS method.
# Stability variable S contributes to A. Element (M, M) of A is filled manually
# as the central diagonal has one more element than the diagonals either side.
A = np.zeros([len(x), len(x)])

for i in range (0, M):
       A[i + 1, i] = S
       A[i, i] = -2 * S + 1
       A[i, i + 1] = S

A[M, M] = -2 * S + 1

# Intialise array for filling with the standard deviation of the diffusing 
# distribution at each time step.
sds = np.zeros(len(t))

# Iterate over each time step to solve FTCS method. Maintain boundary
# conditions from variable uback. Keep track of standard deviation.
for i in range(0, N):
       u[i, 0], u[i, M] = uback, uback
       u[i + 1, :] = np.dot(A, u[i, :])
       sds[i] = sdf(dx, x, u, i, uback)

# Fill in final standard deviation of distribution.
u[N, 0], u[N, M] = uback, uback
sds[N] = sdf(dx, x, u, N, uback)

# Animate the Gaussian being transported.
fig1, ax1 = plt.subplots(figsize = (8, 6))

def animatex(i):
       ax1.clear()
       ax1.set_ylim(uback, np.max(u))
       ax1.plot(x, uback * np.ones(len(x)), label = 'background temperature')
       ax1.plot(x, u[i, :], label = 'total temperature at time t = {:.3f} \
       seconds'.format(round(t[i], decp)))
       ax1.set_xlabel(r'$x$')
       ax1.set_ylabel('temperature')
       ax1.legend()

ani = FuncAnimation(fig1, animatex, repeat = True, interval=50, frames=N + 1)

# Plot time evolution of distribution's standard deviation and compare it to
# analytical solution.
fig2, ax2 = plt.subplots(figsize = (8, 6))

ax2.plot(t, sds, color = 'r', label = 'distribution standard deviation')
ax2.plot(t, np.sqrt(var + (2 * D * t)), color = 'b', label = 'analytical \
standard deviation')
ax2.set_title('Standard deviation 'r'$\sigma$'' as a function of time 'r'$t$')
ax2.set_xlabel(r'$t$'' (s)')
ax2.set_ylabel(r'$\sigma$')
ax2.legend()