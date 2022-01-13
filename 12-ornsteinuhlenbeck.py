"""
Program that performs Ornstien-Uhlenbeck process, using Euler method to solve
first order differential equation.
"""
import random as rd
import numpy as np
import matplotlib.pyplot as plt


"""
Box-Mueller algorithm. Returns a randomly picked number from a Gaussian
distribution, of which the width is a function of variables dt and p.
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
Ornstein-Uhlenbeck Euler method.
"""
def ou(i):
       return x[i] + (l * (mu - x[i]) * dt) + (g * dW(N, sd))


"""
Initialise seed for random number generation, using default argument of time.
"""
rd.seed()

"""
Get user input.
"""
x0 = eval(input("Enter the initial x-position: "))
g = eval(input("Enter the strength of the stochastic term: "))
l = eval(input("Enter the strength of the damping: "))
mu = eval(input("Enter the drift term: "))
tf = eval(input("Enter the final time (in seconds) to iterate to: "))
N = eval(input("Enter the number of intervals to evaluate the function at: "))

"""
Recommended values / values in this region.
x0 = 2
g = 1
l = 1
mu = -2
tf = 20
N = 200
"""

p = 1 / 2
dt = tf / N
sd = np.power(dt, p)

"""
Initialise arrays.
"""
x = np.zeros(N + 1)
x[0] = x0
t = np.zeros(N + 1)
t[0] = 0
Ex = np.zeros(N + 1)
Ex[0] = x[0]
Ex2 = np.zeros(N + 1)
Ex2[0] = np.power(x[0], 2)
Vx = np.zeros(N + 1)
Vx[0] = Ex2[0] - np.square(Ex[0])


"""
Perform Ornstien-Uhlenbeck process.
"""
for i in range(0, N):
       x[i + 1] = ou(i)
       t[i + 1] = (i + 1) * dt
       Ex[i + 1] = sum(x) / (i + 2)
       Ex2[i + 1] = sum(np.power(x, 2)) / (i + 2)
       Vx[i + 1] = Ex2[i + 1] - np.square(Ex[i + 1])


"""
Plotting.
"""
fig, axs = plt.subplots(3, 1, sharex = True)

custom_xlim = (0, t[N])
plt.setp(axs, xlim=custom_xlim)

axs[0].plot(t, x, '+', color = 'g', label = 'Ornstein-Uhlenbeck simulated')
axs[0].plot(t, x0 * np.exp(-l * t) + mu * (1 - np.exp(-l * t)), color = 'r',
            label = 'Noiseless solution')
axs[0].plot(t, x0 * np.exp(-l * t) + mu * (1 - np.exp(-l * t)) + \
            np.sqrt(np.power(g, 2) * (1 - np.exp(-2 * l * t)) / (2 * l)), \
            color = 'y', label = '$\pm$$\sigma$ model')
axs[0].plot(t, x0 * np.exp(-l * t) + mu * (1 - np.exp(-l * t)) - \
            np.sqrt(np.power(g, 2) * (1 - np.exp(-2 * l * t)) / (2 * l)), \
            color = 'y')
axs[0].plot(t, x0 * np.exp(-l * t) + mu * (1 - np.exp(-l * t)) + np.sqrt(Vx), \
            color = 'b', label = '$\pm$$\sigma$ simulated')
axs[0].plot(t, x0 * np.exp(-l * t) + mu * (1 - np.exp(-l * t)) - np.sqrt(Vx), \
            color = 'b')

axs[0].set_title("Time evolution of "r'$x$')
axs[0].set_ylabel(r'$x$')
axs[0].legend()

axs[1].plot(t, Ex, color = 'b', label = r'$\langle x \rangle$')
axs[1].plot(t, mu * np.ones(len(t)), color = 'y', label = 'Long-term mean valu\
e')
axs[1].plot(t, x0 * np.ones(len(t)), color = 'r', label = 'Initial 'r'$x$'' po\
sition')
axs[1].set_title("Time evolution of "r'$\langle x \rangle$')
axs[1].set_ylabel(r'$\langle x \rangle$')
axs[1].legend()

axs[2].plot(t, Vx, color = 'b', label = r'$\sigma^{2}$'' simulated')
axs[2].plot(t, np.power(g, 2) * (1 - np.exp(-2 * l * t)) / (2 * l), \
            color = 'y', label = r'$\sigma^{2}$'' model')
axs[2].set_title('Time evolution of 'r'$\sigma^{2}$')
axs[2].set_ylabel(r'$\sigma^{2}$')
axs[2].set_xlabel(r'$t$')
axs[2].legend()
