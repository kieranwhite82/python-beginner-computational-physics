# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 16:12:41 2021

@author: uqkwhi18
"""

import numpy as np
import matplotlib.pyplot as plt


x = np.linspace(-10, 10, 100)

y = np.sin(x) / x

plt.plot(x, y, 'b*')
plt.xlabel('Time [sec]')
plt.ylabel('Radius [m]')
plt.title('My first plot with matplotlib')
plt.show
