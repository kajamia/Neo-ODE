from Neoode.ode import test_exponential

import numpy as np
import matplotlib.pyplot as plt


all_y = np.asarray(test_exponential())

ie = all_y[:, 0:100]
ee = all_y[:, 100:200]
cn = all_y[:, 200:]

plt.plot(ie[0, :], ie[1, :])
plt.plot(ee[0, :], ie[1, :])
plt.plot(cn[0, :], ie[1, :])

plt.show()
