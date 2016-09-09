import matplotlib.pyplot as plt
import numpy as np

A = np.loadtxt("out.txt")

f1 = A[:,0]
f2 = A[:,1]
h = A[:,2]
#A = np.reshape(A, (np.shape(A)[1], np.shape(A)[0]))
#print(np.shape(A))


plt.semilogx(h, f1, "-")
plt.xlabel("Step size")
plt.ylabel("Error for equation 1")
plt.show()

plt.semilogx(h, f2, "-")
plt.xlabel("Step size")
plt.ylabel("Error for equation 2")
plt.show()