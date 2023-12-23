import subprocess
import numpy as np
import matplotlib.pyplot as plt

# the data will be taken as csv from the output of the command test
test = input("Please enter data source command: ")

input = subprocess.run(test, shell=True, capture_output=True)
print(input.stdout)
data = np.fromstring(input.stdout, sep=",")
t = np.arange(100)

plt.plot(t, data)
plt.show()
