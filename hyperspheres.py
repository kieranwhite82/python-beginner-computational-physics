"""
Program to compute the volume of n-dimensional hyperspheres of three radii.

Input
Radius of hyperspheres.
Maximum dimension nmax of the hyperspheres (integer).

Output
Volume of each hypersphere.
Plot the hypersphere's volume as a function of n.
CSV file of data.
"""
import numpy as np
import matplotlib.pyplot as plt
import csv
from pandas import read_csv


# Function to calculate single factorial.
def nf(n):
    i = 0
    ni = n
    if n > 1:
        while ni - i > 0:
            n *= ni - i
            i += 1
        return n/ni
    elif n == -1 or n == 0 or n == 1:
        return 1.0
    else:
        print("The factorial function failed.")


# Function to calculate double factorial.
def ndf(n):
    i = 0
    ni = n
    if n > 0 and n % 2 == 0:
        while ni - (2*i) > 1:
            n *= ni - (2*i)
            i += 1
        return n/ni
    elif n > 0 and n % 2 == 1:
        while ni - (2*i) > 0:
            n *= ni - (2*i)
            i += 1
        return n/ni
    elif n == -1 or n == 0:
        return 1.0
    else:
        print("The double factorial function failed.")


# Function to calculate the volume of a hypersphere.
def V(n, r):
    x = np.zeros(n+1)
    y = np.zeros(n+1)
    while n >= 0:
        if n % 2 == 0 and n > 0:
            y[n] = np.power(2, 1) * np.power(np.pi, n/2) * np.power(r, n) /\
                (n * nf((n/2)-1))
            x[n] = n
            n -= 1
        elif n % 2 == 1 and n > 0:
            y[n] = np.power(2, (n+1)/2) * np.power(np.pi, (n-1)/2) *\
                np.power(r, n) / (n * ndf(n-2))
            x[n] = n
            n -= 1
        elif n == 0:
            x[0] = 0
            y[0] = 1
            n -= 1
        else:
            print("Invalid maximum dimension entered.")
    return x, y


# User input for variables.
nmax = eval(input("Enter the maximum dimension of the hyperspheres: "))
radius1 = eval(input("Enter a first radius of the hypersphere: "))
radius2 = eval(input("Enter a second radius of the hypersphere: "))
radius3 = eval(input("Enter a third radius of the hypersphere: "))
print("\nnmax! = {}".format(nf(nmax)))
print("\nnmax!! = {}".format(ndf(nmax)))


# Calculate volumes for each radius.
ndim1, nvol1 = V(nmax, radius1)
ndim2, nvol2 = V(nmax, radius2)
ndim3, nvol3 = V(nmax, radius3)


# Create matrix to store relevant values for writing to csv data file.
datarray = np.matrix.transpose(np.array([ndim1, nvol1, nvol2, nvol3]))


# Open, write to and close csv data file.
filename = "6-data.csv"
f = open(filename, "w", newline='')
wr = csv.writer(f)
wr.writerows(datarray)
f.close()


# Add headers to columns.
df = read_csv('6-data.csv', header=None)
df.columns = ['n', 'Vn(r={})'.format(radius1), 'Vn(r={})'.format(radius2), 'Vn(r={})'.format(radius3)]  # noqa
df.to_csv('6-data.csv', index=False)


# Plot resulting hypersphere volumes from user input.
plt.figure()
plt.plot(ndim1, nvol1, 'bx', label='Radius {}'.format(radius1))
plt.plot(ndim2, nvol2, 'rx', label='Radius {}'.format(radius2))
plt.plot(ndim3, nvol3, 'yx', label='Radius {}'.format(radius3))
plt.xlabel("n")
plt.ylabel("Vn")
plt.legend()
plt.title("Hypersphere volume up to {} dimension, of radii {}, {}, {}.".format(nmax, radius1, radius2, radius3))  # noqa
