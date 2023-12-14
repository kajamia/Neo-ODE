from Neoode.ode import test_mass_spring, test_alpha
from Neosoft.cla import Matrix
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



# generic alpha

all_y = test_alpha()

r1 = all_y[:, 0]
r2 = all_y[:, 1]
r3 = all_y[:, 2]

plt.plot(t, r1)
plt.plot(t, r2)
plt.plot(t, r3)

plt.show()
