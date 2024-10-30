import numpy as np

import matplotlib.pyplot as plt

# Read the zmForward matrix from file
zmForward = np.loadtxt('./zmForward.txt')

# Plot the matrix using plt.imshow
plt.imshow(zmForward, cmap='viridis')
plt.colorbar()
plt.show()