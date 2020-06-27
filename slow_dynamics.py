from Network import Network
from helpers_slow import power_vector
import matplotlib.pyplot as plt
import numpy as np

# global vars
NR_GENERATORS = 12
ABS_FLOW = []
MAX_FLOW = []
MAX_OVERLOAD = []

def initialize_model(n_buses, max_P_gen, layers, gen_layer, max_lambda, days):
	'''
	PARAMETERS:
	n_buses: number of buses
	max_P_gen:
	layers: Total number of layers in the network
	gen_layer: The layer at which the generators will be placed
	max_lambda: The maximum value of lambda
	days: The number of days for the simulation

	RETURNS:
	F_0: The flow on day 0
	A: matrix which represents the network constraints
	P_0: The power injection on day 0
	daily_lambda: daily multiplication, it represents a slowly increasing secular load
	max_F: The limit of flow

	'''
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

def simulation(days, A, P_0, daily_lambda, max_F, F_0):
	'''
	PARAMETERS:
	days:  The number of days for the simulation
	A: matrix which represents the network constraints
	P_0: The power injection on day 0
	daily_lambda: daily multiplication, it represents a slowly increasing secular load
	F_0: The flow on day 0
	max_F: The limit of flow


	'''

	mean_P_load = []
	mean_P_gen = []
	max_flow_gen = []

	P_k = P_0

	for k in range(days):

		# increase the power in the buses with lambda and update flow
		P_k = P_k * daily_lambda[k]
		F_k = A.dot(P_k)


if __name__ == '__main__':
	n_buses = 94
	max_P_gen = 10
	layers = 6
	gen_layer = 3

	lambda_increase = 0.0005
	mean_lambda = 1 + lambda_increase
	max_lambda = 1 + 2 * lambda_increase

	days = 3
	F_0, A, P_0, daily_lambda, max_F = initialize_model(n_buses, max_P_gen, layers, gen_layer, max_lambda, days)
    simulation(days, A, P_0, daily_lambda, max_F, F_0)
