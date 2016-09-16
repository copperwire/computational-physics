import matplotlib.pyplot as plt
import numpy as np
import sys

for i in range(5):
	try: 
		A_gen = np.loadtxt("forw_back" + "{}".format(i))
		A_spec = np.loadtxt("special_result"+ "{}".format(i))
		
		x_spec = A_spec[:,0]
		U_spec = A_spec[:,1]
		exact_spec = A_spec[:,2]
		eps_spec = A_spec[:,3]

		x_gen = A_spec[:,0]
		U_gen = A_spec[:,1]
		exact_gen = A_spec[:,2]
		eps_gen = A_spec[:,3]

		Np1 = len(x_gen)
		
		plt.plot(x_gen, exact_gen, "r-", label = "Exact solution")
		plt.plot(x_gen, U_gen, "b-", label = "gen algor.")
		plt.plot(x_spec, U_spec, "g-", label = "spec algor.")
		plt.legend()
		plt.title("Exact vs numerical algorithms")
		plt.ylabel(r"$y_i$")
		plt.xlabel(r"$x_i$")
		plt.savefig("figure_"+"{}".format(i))
		plt.cla()
		plt.clf()

		plt.plot(x_gen[1:Np1-1], eps_gen[1:Np1-1], label = "General algorithm")
		plt.plot(x_spec[1:Np1-1], eps_spec[1:Np1-1], label = "Special algorithm")
		plt.legend()
		plt.title(r"Relative error, $\epsilon_r$, comparison")
		plt.ylabel(r"$\epsilon_r$")
		plt.xlabel(r"$x_i$")
		plt.savefig("error_plot_{}".format(i))
		plt.cla()
		plt.clf()

		"latex writing code should go here"
	except FileNotFoundError:
		continue