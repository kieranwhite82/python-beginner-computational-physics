"""
Plot sin function of radial position as a function of time.
(Circular motion)
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
