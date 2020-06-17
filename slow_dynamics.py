from make_network_and_matrices import Network
from helpers_slow import power_vector
import matplotlib.pyplot as plt
import numpy as np

# global vars
NR_GENERATORS = 12

# variables
n_buses = 94
max_P_gen = 10
max_F = 6

# total nr of layers in the network
layers = 6

# the layer at which the generators will be placed
gen_layer = 3


def initialize_model(n_buses, max_P_gen, layers, gen_layer, max_F):

	# create the network and matrix A
	N = Network(layers, gen_layer)
	# print(N.A)
	N.plot_network()
	A = N.matrix_a()

	# create the power vector at day 0
	P_0 = power_vector(n_buses, max_P_gen, NR_GENERATORS)

	if abs(sum(P_0)) > 0.0001:
		print("SUM P not approx 0 -> ending initialization")
		exit()

	for value in P_0:
		print("val", value)
	# calculate the flow on each line at day 0 based on matrix A and the power vector
	F_0 = A.dot(P_0)

	for flow in F_0:
		print(flow)
		# print(len(flow))
		if abs(flow) > max_F:
			print("Flow at day 0 too large -> ending initialization")
			exit()

	return A, P_0


A, P_0 = initialize_model(n_buses, max_P_gen, layers, gen_layer, max_F)
