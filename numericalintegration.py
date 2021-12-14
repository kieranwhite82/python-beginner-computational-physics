"""
Program that applies three numerical integration methods to solve integrand
from fixed list of possibilities. 
Updated to also take user input function with appropriate numpy functions.
Program can calculate net area or absolute area.

Input
Net or absolute area.
Integrand from fixed list of possiblities or user input function.
Integration limits.
Desired accuracy.

Output
Estimate for integral once the user input accuracy has been reached for each
method.
"""
import numpy as np


# Exit functions for invalid selections.
def f1(): raise Exception("Invalid number/character entered for integrand selection.")

def f2(): raise Exception("Invalid number/character entered for area calculation type.")

def f3(): raise Exception("The user entered limits were numerically equal, therefore the integral is equal to zero.")
        

# Get user input for the integrand.
print("Integrand option list:\nFor (1+x^2)^-1, enter 1\nFor 1-x^2, enter 2\nFor 1/x, enter 3\nTo enter user defined function, enter 4")
Is = eval(input("For net area, enter 1. For absolute area, enter 2: "))
fI = eval(input("Enter the number corresponding to the integrand from the option list to numerically integrate: "))  # noqa
if fI > 0 and fI < 4 and (fI % 2 == 0 or fI % 2 == 1):
    print("Integrand {} selected.".format(fI))
elif fI == 4:
    userequ = input("Enter a function of x (using numpy mathematical functions, using numpy as np, for example: np.exp(x); the same applies to special limits, for example: np.pi): ")
else:
    f1()


# Assign relevant integrand to function to numerically integrate.
if Is == 1:
    if fI == 1:
        def y(x):
            return np.power(1 + np.power(x, 2, dtype=float), -1, dtype=float)
    elif fI == 2:
        def y(x):
            return 1 - np.power(x, 2, dtype=float)
    elif fI == 3:
        def y(x):
            return np.power(x, -1, dtype=float)
    else:
        def y(x):
            return eval(userequ)
elif Is == 2:
    if fI == 1:
        def y(x):
            return abs(np.power(1 + np.power(x, 2, dtype=float), -1, dtype=float))
    elif fI == 2:
        def y(x):
            return abs(1 - np.power(x, 2, dtype=float))
    elif fI == 3:
        def y(x):
            return abs(np.power(x, -1, dtype=float))
    else:
        def y(x):
            return abs(eval(userequ))
else:
     f2()


# Define numerical integration functions to call on in numerical solver loops.
# Rectangular rule:
def RR(a, h, N):
    f = 0
    F = 0
    i = 0
    while i < N:
        f = h * y(a + ((i + 0.5) * h))
        F += f
        i += 1
    return F


# Trapezoidal rule:
def TR(a, b, h, N):
    f = 0
    F = h / 2 * (y(a) + y(b))
    i = 1
    while i < N:
        f = h * (y(a + (i * h)))
        F += f
        i += 1
    return F


# Simpson's rule:
def SR(a, b, h, N):
    f = 0
    F = h / 3 * (y(a) + y(b))
    i = 1
    while i < N:
        if i % 2 == 1:
            f = 4 * h / 3 * (y(a + (i * h)))
            F += f
            i += 1
        elif i % 2 == 0:
            f = 2 * h / 3 * (y(a + (i * h)))
            F += f
            i += 1
        else:
            print("Even/odd summation failed.")
    return F


# Get user input for limits to integral.
a = eval(input("Enter the lower limit to the integral: "))
b = eval(input("Enter the upper limit to the integral: "))

# Dummy variables for switching limits in case of user entering b < a.
c = 0
d = 0

# Response to user input limits.
if b == a:
    f3()
elif a > b:
    c = a
    d = b
    a = d
    b = c
    print("User entered limits are reversed in the case of the upper limit being numerically smaller than the lower limit.")
else:
    print("User entered suitable integral limits.")


# The number of intervals N starts at 10, and will double until the user
# entered accuracy E is reached. The interval spacing h is calculated given
# user input.
N = 10
h = (b - a) / N

# The accuracy here is defined as the difference between successive integration
# results, the idea being as the iterations approaches a value, the difference
# between the results will decrease.
E = eval(input("Enter the desired accuracy for the integral: "))


# Rectangular rule iteration.
print("\nRectangular rule numerical integration:")
FRRn = 0
NRR = N
hRR = h
ERR = E + 1
while E < ERR:
    FRRo = FRRn
    FRRn = RR(a, hRR, NRR)
    ERR = abs(FRRn - FRRo)
    NRR = 2 * NRR
    hRR = (b - a) / NRR

print("The accuracy {} to iterate to was reached with the numerical integral as {} with {} integration steps.".format(E, FRRn, NRR))


# Trapezoidal rule iteration.
print("\nTrapezoidal rule numerical integration:")
FTRn = 0
NTR = N
hTR = h
ETR = E + 1
while E < ETR:
    FTRo = FTRn
    FTRn = TR(a, b, hTR, NTR)
    ETR = abs(FTRn - FTRo)
    NTR = 2 * NTR
    hTR = (b - a) / NTR

print("The accuracy {} to iterate to was reached with the numerical integral as {} with {} integration steps.".format(E, FTRn, NTR))


# Simpson's rule iteration.
print("\nSimpson's rule numerical integration:")
FSRn = 0
NSR = N
hSR = h
ESR = E + 1
while E < ESR:
    FSRo = FSRn
    FSRn = SR(a, b, hSR, NSR)
    ESR = abs(FSRn - FSRo)
    NSR = 2 * NSR
    hSR = (b - a) / NSR

print("The accuracy {} to iterate to was reached with the numerical integral as {} with {} integration steps.".format(E, FTRn, NTR))
