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

# #Define M for day k
#TODO: loop over days to get M
days = 10
for k in days:
    A, P_0 = initialize_model(n_buses, max_P_gen, layers, gen_layer, max_F)
    Fk = np.dot(A,P_0)
    Mk = [i/max_F for i in Fk]
    P_outaged = [0.001*i for i in Mjk]
    P_outaged_sum = sum(P_outaged)
    P_outaged = [i/P_outaged_sum for i in P_outaged] #give probabilities
    P_overload = [0.3 *i for i in Mjk]
    P_overload_sum = sum(P_outaged)
    P_overload = [i/P_overload_sum for i in P_overload] #give probabilities

    print("P_outaged is {},\n P_overload is {},\n".format(P_outaged,P_overload))
