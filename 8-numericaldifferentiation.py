# -*- coding: utf-8 -*-
"""
Created on Sun Dec 12 20:03:25 2021

@author: kiera

Program that performs numerical differentiation: three methods for first order
ODE's, and one method for second order ODE's.

Input
User input function with appropriate numpy functions.
Initial value for IVP x0 and end point xe.

Output
Numerical differentiation results from the three first order ODE methods, and
one result from the second order ODE method.
"""
import numpy as np
import matplotlib.pyplot as plt

print("\nThis program attempts to numerically differentiate user-input \
functions. It helps when functions are smooth over the interval: functions \
such as tangent function near n*pi/2 where n is an integer approach infinity, \
and so are not appropriate for this program. Another example is dirac delta \
function.")

print("\nThis program uses numpy mathematical functions, so enter functions \
to numerically differentiate appropriately.")

# Get user input for the function to numerically differentiate.
userequ = input("Enter a function of x: ")


def f(x):
       return eval(userequ)


# Get user input for bounds within which to numerically differentiate.
x0 = eval(input("Enter the initial value at which to evaluate the function: "))
xf = eval(input("Enter the final value at which to evaluate the function: "))

# Get user input on the number of intervals at which to estimate the
# differential function.
N = eval(input("Enter the number of intervals at which to estimate the \
differential function: "))

# Calculate step size:
dx = (xf - x0) / N

# Setup arrays for storing estimates for differential function at each step:
x = np.zeros(N + 1)
F = np.zeros(N + 1)


# Forward finite difference function:
def ffd(x0, i, dx):
       return (f(x0 + ((i + 1) * dx)) - f(x0 + (i * dx))) / dx


# Central finite difference function:
def cfd(x0, i, dx):
       return (f(x0 + ((i + 1) * dx)) - f(x0 + ((i - 1)* dx))) / (2 * dx)


# 4-point central difference approximation function:
def pcd(x0, i, dx):
       return ( - f(x0 + ((i + 2) * dx)) + f(x0 + ((i - 2)* dx)) + (8 * (f(x0 \
               + ((i + 1) * dx)) - f(x0 + (i * dx))))) / (4 * dx)

# Correction factor of 3 applied to pcd function. Not sure why values are
# initally scaled by 1/3.


# Central finite difference, second order function:
def cfd2(x0, i, dx):
       return (f(x0 + ((i + 1) * dx)) - (2 * f(x0 + (i * dx))) + f(x0 + ((i - \
               1)* dx))) / (dx * dx)


# Initialise arrays, variables for looping
i = 0

dFffd = np.zeros(N + 1)
dFcfd = np.zeros(N + 1)
dFpcd = np.zeros(N + 1)
ddFcfd = np.zeros(N + 1)


# Calculate approximations for derivatives:
while i < N + 1:
       x[i] = x0 + (i * dx)
       F[i] = f(x0 + (i * dx))
       dFffd[i] = ffd(x0, i, dx)
       dFcfd[i] = cfd(x0, i, dx)
       dFpcd[i] = pcd(x0, i, dx)
       ddFcfd[i] = cfd2(x0, i, dx)
       i += 1


# Plot derivative approximations.
plt.figure()
plt.plot(x, F, 'b', label='F(x)')
plt.plot(x, dFffd, 'r', label='ffd')
plt.plot(x, dFcfd, 'g', label='cfd')
plt.plot(x, dFpcd, 'c', label='pcd')
plt.plot(x, ddFcfd, 'y', label='cfd2')
plt.legend()
