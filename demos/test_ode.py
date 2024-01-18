from Neoode.ode import *
from Neosoft.cla import Matrix, Vector
import numpy as np
import matplotlib.pyplot as plt


# explicit and implicit euler

all_y = np.asarray(test_mass_spring())

t = np.linspace(0, 4*np.pi, 100)

ie = all_y[0:100, :]
ee = all_y[100:200, :]
cn = all_y[200:, :]

plt.plot(t, ie[:, 1], "b")
plt.plot(t, ie[:, 1], "g")
plt.plot(t, ie[:, 1], "r")

plt.plot(t, ie[:, 0], "y")
plt.plot(t, ee[:, 0], "orange")
plt.plot(t, cn[:, 0], "black")

plt.show()
