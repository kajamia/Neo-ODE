#! /bin/python

import subprocess
import sys
import numpy as np
import matplotlib.pyplot as plt

# the data will be taken as csv from the output of the command passed as first argument
test = sys.argv[1]

input = subprocess.run(test, shell=True, capture_output=True)
print(input.stdout)
data = np.fromstring(input.stdout, sep=",")
t = np.arange(len(data))

plt.plot(t, data)
plt.show()
