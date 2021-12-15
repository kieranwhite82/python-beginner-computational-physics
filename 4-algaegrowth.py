"""
Program that computes the size of an algae bloom after some time t.

Input
Initial bloom size.
Bloom size units.
Growth constant.
Growth constant units.
Time t in the future to calculate bloom size.
The size of the lake the bloom is on.

Output
Predicted bloom size at the time t.
The time when the lake will be completely covered by the bloom.
"""
import numpy as np


s0 = eval(input("Enter the intial bloom size value: "))
sU = input("Enter the length units used for the initial bloom size: ")
print("\nThe initial bloom size is {} {}^2.".format(s0, sU))

a = eval(input("Enter the growth constant value: "))
tU = input("Enter the units of time used for bloom growth calculation (\
inverse of units of growth constant): ")
print("\nThe growth constant is {} {}^-1.".format(a, tU))

lA = eval(input("Enter the size of the lake the bloom is on with length units \
{}: ".format(sU)))
tl = np.log(lA / s0) / a

t = eval(input("Enter the time at which to calculate the size of the bloom in \
{}: ".format(tU)))
if t < tl:
    s = s0 * np.exp(a * t)
    print("\nThe size of the bloom after {} {} is {} {}^2. It will take \
    ~{:.1f} {} for the bloom to cover the lake.".format(t, tU, s, sU, tl, tU))
else:
    print("\nThe time entered to calculate the bloom size at was greater than \
    the time taken for the bloom's area to grow to the lake's area. It will \
    take ~{:.1f} {} for the bloom to cover the lake.".format(tl, tU))
