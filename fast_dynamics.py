from make_network_and_matrices import Network
from helpers_slow import power_vector
import matplotlib.pyplot as plt
import numpy as np
from slow_dynamics import initialize_model

# global vars
NR_GENERATORS = 12
# variables
n_buses,max_P_gen,max_F= 94,10,6
# total nr of layers in the network
layers = 6
# the layer at which the generators will be placed
gen_layer = 3

A, P_0 = initialize_model(n_buses, max_P_gen, layers, gen_layer, max_F)
Fk = np.dot(A,P_0)
# Fjk = np.zeros([Fk])
# Mjk = []

# print(A.shape)
# #Define M
# for i in range(len(F)):
Mjk = [i/max_F for i in Fk]
P_outaged = [3*i + 5 for i in Mjk]
P_outaged_sum = sum(P_outaged)
P_outaged = [i/P_outaged_sum for i in P_outaged] #give probabilities
print(P_outaged)
# P_overload= [3*i + 5 for i in Mjk]

# print(P_outaged)
