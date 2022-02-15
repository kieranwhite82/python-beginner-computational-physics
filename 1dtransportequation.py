"""
Program that solves the transport PDE through the FTCS method with the initial
condition of a gaussian. The boundary conditions are zero. Lax-Friedrich scheme
is applied to the FTCS method. One-dimensional space.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Velocity v is the constant velocity at which the transported function moves
# along the position axis.
v = 8

# Time t: t0 =< t =< tf. Have N steps at which to evaluate the FTCS scheme. The
# time interval is dt.
t0 = 0
tf = 2.5
N = 100
dt = tf / N
t = np.linspace(t0, tf, num = N + 1, endpoint = True)

# Position x: x0 =< x =< xf. There are M total discrete positions. There are
# the boundary conditions at the ends: x0  and xf. Therefore M - 2 positions
# are to be evaluated, with position interval dx.
x0 = 0
xf = 20
M = 100
dx = xf / M
x = np.linspace(x0, xf, num = M + 1, endpoint = True)

# The gaussian function u(x, t) is transported along the position axis. It has
# normalisation constant Nc, and width var. The inital position of the
# function is u0 = u(x, t0). The inital offset with magnitude 4 * sd is there
# to conserve smoothness requirements for the boundary conditions.
var = 1
sd = np.sqrt(var)
nc = 1 / np.sqrt(2 * np.pi * var)
u0 = nc * np.exp(-np.power(x - (4 * sd), 2) / (2 * var))
u = np.zeros([len(t), len(x)])
u[0, :] = u0
u[:, 0], u[:, M] = 0, 0

# For reporting on stability condition, S. If u(x, t0), u(x0, t) and u(xf, t)
# are smooth, and S =< 1, the method is stable.
S = np.abs(v * dt / dx)

print("The stability variable is {}".format(S))
if S == 1:
       print("The stability variable is equal to unity, shape should be \
conserved to boundaries and solution should be stable.")
elif S < 1:
       print("The stability variable is less than unity, shape will not be \
conserved to boundaries and solution will be stable.")
else:
       print("The stability variable is greater than unity, shape will not be \
conserved to boundaries and solution will not be stable.")

# Tridiagonal matrix A to solve set of linear equations making up FTCS method.
# Simplification variable a contributes to A.
A = np.zeros([len(x) - 2, len(x) - 2])
a = (v * dt) / (2 * dx)
for i in range (0, M - 2):
       A[i + 1, i] = 0.5 + a
       A[i, i + 1] = 0.5 - a

# Iterate over each time step to solve FTCS method.
for i in range(0, N):
       u[i + 1, 1:M] = np.dot(A, u[i, 1:M])

# Animate the Gaussian being transported.
fig1, ax1 = plt.subplots(figsize=(8, 6))

def animatex(i):
       ax1.clear()
       ax1.set_ylim(0, np.max(u))
       ax1.plot(x, u[i, :])

ani = FuncAnimation(fig1, animatex, repeat = True, interval=50, frames=N + 1)
