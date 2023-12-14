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

all_y = np.array(100, 3)
tend = 2*2*np.pi
steps = 100
x  = Vector(3)
dx  = Vector(3)
ddx  = Vector(3)

rhs = dLagrange()
mass = Projector(3, 0, 2)

for i in range(100):
    SolveODE_Alpha(tend, steps, 0.8, x, dx, ddx, rhs, mass)
    
