"""
Program that generates a pseudo-random sample with a uniform distribution, and
transforms it using Box-Mueller algorithm to generate a Gaussian distribution, 
one-dimensional and two-dimensional.

Input
Number of inital pseudorandom sample points to simulate
Expected value of Gaussian distribution (mean)
Standard deviation
Number of standard deviations to bin sample points into
Decimal places to round standard deviation to for representation

Output
Uniform one-dimensional (bar) and two-dimensional (scatter) distributions
Gausian one-dimensional (bar) and two-dimensional (scatter) distributions
"""
import random as rd
import numpy as np
import matplotlib.pyplot as plt

# Initialise seed fo rrandom number generation, using default argument of time.
rd.seed()

# Get user input on variables.
N = eval(input("Enter the number of inital sample points: "))
mu1 = eval(input("Enter the expectation value of the Gaussian distribution: "))
sd = eval(input("Enter the standard deviation: "))
nsd = eval(input("Enter the number of standard deviations from the mean to bin\
 sample points into: "))
sdR = eval(input("Enter the number of decimal points to round sigma to for fig\
ures: "))

# Construct lists of bin spacings. For the initial sample, it is split into
# fourths. For the Gaussian-transformed points, the spacing is the standard
# deviation.
a = []
for i in range(0, 4):
       
       A = i / 4
       a.append(A)

x = []
for i in range(0, (2 * nsd)):
       
       X = sd * (-nsd + i)
       x.append(X)


# Box-Mueller algorithm.
ic = 0
r1 = np.zeros(N)
r2 = np.zeros(N)
r1s = np.zeros(N)
r2s = np.zeros(N)
r = np.zeros(N)
j = 0
s = []
x1 = []
x2 = []
for i in range(0, N):
       
       r1[i] = rd.random()
       r2[i] = rd.random()
       r1s[i] = (2 * r1[i]) - 1.0
       r2s[i] = (2 * r2[i]) - 1.0
       r[i] = r1s[i]**2 + r2s[i]**2

       if r[i] > 0 and r[i] < 1:
              S = np.sqrt(-2 * np.log(r[i]) / r[i])
              s.append(S)
              X1 = (r1s[i] * s[j] * sd) + mu1
              x1.append(X1)
              X2 = r2s[i] * s[j] * sd
              x2.append(X2)
              j += 1
       else:
              pass


# Return the number of initial sample points that were appropriate during
# transformation process to be Gaussian-distributed.
print("\n{} of the initial {} sample points were valid for Gaussian transforma\
tion.".format(len(s), N))

# Bin points for both sample groups.
a1 = np.zeros(4)
for i in range (0, N):
       
       for n in range(0, 4):
              
              if r1[i] >= n * 1 / 4 and r1[i] < (1 / 4) * (n + 1):
                     a1[n] += 1
              else:
                     pass

b1 = np.zeros(2 * nsd)
for i in range (0, len(s)):
       
       for n in range(0, 2 * nsd):
              
              if x1[i] - mu1 >= ((-nsd + n) * sd) and x1[i] - mu1 < ((-nsd + n\
                                                      + 1) * sd):
                     b1[n] += 1
              else:
                     pass
Ex = sum(x1) / len(s)
Ex2 = sum(np.power(x1, 2)) / len(s)
sdc = np.sqrt(Ex2 - np.power(Ex, 2))
print("\n{} of the total {} Gaussian distributed points were within +-{} stand\
ard deviations.".format(int(sum(b1)), len(s), round(nsd * sd, sdR)))
print("\nThe expectation value of the Gaussian-distributed variables is {}; it\
 should approach {}. The expectation value of the square of the Gaussian-distr\
ibuted values is {}. This results in a calculated standard deviation of {}, co\
mpared to the user entered value {}.".format(Ex, mu1, Ex2, sdc, sd))

# Plot results.
fig, axs = plt.subplots(2, 2)

axs[0, 0].bar(a, a1 / N, width = 1 / 4, align = 'edge')
axs[0, 0].set_xticks(np.arange(0, 1 + 1 / 4, 1 / 4))
axs[0, 0].set_xticklabels(np.arange(0, 1 + 1 / 4, 1 / 4).tolist())
axs[0, 0].set_title('Discrete uniform probability distribution function')
axs[0, 0].set_xlabel('Pseudorandom number E (0, 1)')
axs[0, 0].set_ylabel('Fraction')

axs[0, 1].plot(r1, r2, '+')
axs[0, 1].set_xticks(np.arange(0, 1 + 1 / 4, 1 / 4))
axs[0, 1].set_xticklabels(np.arange(0, 1 + 1 / 4, 1 / 4).tolist())
axs[0, 1].set_yticks(np.arange(0, 1 + 1 / 4, 1 / 4))
axs[0, 1].set_yticklabels(np.arange(0, 1 + 1 / 4, 1 / 4).tolist())
axs[0, 1].set_title('2-D uniform sample point distribution')
axs[0, 1].set_xlabel('Pseudorandom number X E (0, 1)')
axs[0, 1].set_ylabel('Pseudorandom number Y E (0, 1)')

axs[1, 0].set_title('Discrete Gaussian probability distribution function')
axs[1, 0].set_xlabel('Gaussian variable X')
axs[1, 0].set_ylabel('Fraction')

lbl = []
for i in range(0, (2 * nsd) + 1):
       if i == nsd:
              LBL = r"$\mu$"
              lbl.append(LBL)
       else:
              LBL = r"${0}\sigma$".format(round(sd * (-nsd + i), sdR))
              lbl.append(LBL)

s_ticks = np.arange(-nsd * sd, (nsd * sd) + sd, sd)

axs[1, 0].set_xticks(s_ticks)
axs[1, 0].set_xticklabels(lbl)
axs[1, 0].bar(x, b1 / len(s), width = sd, align = 'edge')
axs[1, 0].set_xlim(-nsd * sd, nsd * sd)

axs[1, 1].plot(x1 - (mu1 * np.ones(len(s))), x2, '+')
axs[1, 1].set_title('2-D Gaussian sample distribution')
axs[1, 1].set_xlabel('Gaussian variable X')
axs[1, 1].set_ylabel('Gaussian variable Y')
axs[1, 1].set_xlim(-nsd * sd, nsd * sd)
axs[1, 1].set_xticks(s_ticks)
axs[1, 1].set_xticklabels(lbl)
axs[1, 1].set_ylim(-nsd * sd, nsd * sd)
axs[1, 1].set_yticks(s_ticks)
axs[1, 1].set_yticklabels(lbl)
