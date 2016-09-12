import matplotlib.pyplot as plt
import numpy as np

U = np.loadtxt("forw_back.txt")
U = U[::-1]
N = np.size(U)

exact = lambda x: 1-(1-np.exp(-10))*x - np.exp(-10*x)

x = np.linspace(0.001, 0.999, N - 2)
X = np.linspace(0, 1, N)
u = np.zeros(N)

u[1:N-1] = exact(x)


plt.plot(x, u[1:N-1], "r-", label = "Exact solution")
plt.plot(X, U, "b-x", label = "Forw + back subst")
plt.legend()
plt.show()