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



# generic alpha

""" all_y = np.array(100, 3)
tend = 2*2*np.pi
steps = 100
x  = Vector(3)
dx  = Vector(3)
ddx  = Vector(3)

rhs = dLagrange()
mass = Projector(3, 0, 2)

for i in range(100):
    SolveODE_Alpha(tend, steps, 0.8, x, dx, ddx, rhs, mass) """

tend = 2*2*np.pi
steps = 100

x = [Vector(3) for _ in range(100)]
dx = [Vector(3) for _ in range(100)]
ddx = [Vector(3) for _ in range(100)]

rhs = dLagrange()
mass = Projector(3, 0, 2)

for i in range(100):
    x[i][0] = 1
    x[i][1] = 0
    x[i][2] = 0

    dx[i][0] = 0
    dx[i][1] = 0
    dx[i][2] = 0

    ddx[i][0] = 0
    ddx[i][1] = 0
    ddx[i][2] = 0

    SolveODE_Alpha(tend*i/100, steps, 0.8, x[i], dx[i], ddx[i], rhs, mass)


y = np.array([y(0) for y in x])
t = np.arange(100)

plt.plot(t, y)
