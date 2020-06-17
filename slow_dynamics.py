from Network import Network
from helpers_slow import power_vector
import matplotlib.pyplot as plt
import numpy as np

# global vars
NR_GENERATORS = 12
ABS_FLOW = []
MAX_FLOW = []
MAX_OVERLOAD = []

# variables
n_buses = 94
max_P_gen = 10

# total nr of layers in the network
layers = 6

# the layer at which the generators will be placed
gen_layer = 3

# average power increase per day
lambda_increase = 0.0005
mean_lambda = 1 + lambda_increase
max_lambda = 1 + 2 * lambda_increase

# the nr of days for the simulation
days = 1000


def initialize_model(n_buses, max_P_gen, layers, gen_layer, max_lambda, days):

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

	# calculate the flow on each line at day 0 based on matrix A and the power vector
	F_0 = A.dot(P_0)

	# define the maximum flow for each line
	max_F = [np.max([abs(flow) for flow in F_0]) * 1.05] * n_buses

	daily_lambda = (max_lambda - 1) * np.random.random_sample(days) + 1


	return F_0, A, P_0, daily_lambda, max_F


F_0, A, P_0, daily_lambda, max_F = initialize_model(n_buses, max_P_gen, layers, gen_layer, max_lambda, days)

def simulation(days, A, P_0, daily_lambda, max_F, F_0):

	mean_P_load = []
	mean_P_gen = []
	max_flow_gen = []


	P_k = P_0


	for k in range(days):

		# increase the power in the buses with lambda and update flow
		P_k = P_k * daily_lambda[k]
		F_k = A.dot(P_k)

		# below this is all to see if the parameters are good
		# P_load = []
		# P_gen = []
		# for value in P_k:
		# 	if value <= 0:
		# 		P_load.append(value)
		# 	else:
		# 		P_gen.append(value)
		# mean_P_load.append(np.mean(P_load))
		# mean_P_gen.append(np.mean(P_gen))
		# ABS_FLOW.append(np.mean([abs(flow) for flow in F_k]))
		# MAX_FLOW.append(np.max([abs(flow) for flow in F_k]))
		# M_k = [abs(flow/max_F[i]) for i, flow in enumerate(F_k)]
		# MAX_OVERLOAD.append(np.max(M_k))



	# plot the flow and overload values to tune params
	# plt.plot(range(0, len(ABS_FLOW)), ABS_FLOW, label = "mean")
	# plt.plot(range(0, len(MAX_FLOW)), MAX_FLOW, label = "max")
	# plt.legend()
	# plt.show()
	# plt.plot(range(0, len(MAX_OVERLOAD)), MAX_OVERLOAD, label = "max overload")
	# plt.legend()
	# plt.show()

simulation(days, A, P_0, daily_lambda, max_F, F_0)

