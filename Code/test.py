import numpy as np

l = [0, 1, 2, 3, 5]

z = np.zeros(7)

print(z)

#z[l] = 1

np.put(z, l, 1)

print(z)

