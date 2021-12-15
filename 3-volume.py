# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 16:29:23 2021

@author: uqkwhi18

Program to calculate the volume of sphere.

Input
Radius of sphere.

Output
Volume of sphere.
"""
import numpy as np


rsphere = eval(input("Enter the radius of the sphere in metres: "))

volume = 4 * np.pi * np.power(rsphere, 3) / 3

print("The volume of a sphere with radius", rsphere, "m is", '%.2f' % volume, "m3.")
