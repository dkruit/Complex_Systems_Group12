from Network import Network
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

# average power increase per day
lambda_increase = 0.0005
mean_lambda = 1 + lambda_increase
max_lambda = 1 + 2 * lambda_increase

# #Define M for day k
#TODO: loop over days to get M
days = 1
# for k in days:
F_0, A, P_0, daily_lambda, max_F = initialize_model(n_buses, max_P_gen, layers, gen_layer, max_lambda, days)

def simulation(days, A, P_0, daily_lambda, max_F, F_0):

    mean_P_load = []
    mean_P_gen = []
    max_flow_gen = []


    P_k = P_0
    # F_max = max_F

    for k in range(days):
        # increase the power in the buses with lambda and update flow
        P_k = P_k * daily_lambda[k]
        F_k = A.dot(P_k)
        F_max = max(F_k)
        M_k = [i/F_max for i in F_k]
        P_outaged = [0.001*i for i in M_k]
        P_outaged_sum = sum(P_outaged)
        P_outaged = [abs(i/P_outaged_sum) for i in P_outaged] #give probabilities
        P_overload = [0.3 *i for i in M_k]
        P_overload_sum = sum(P_outaged)
        P_overload = [abs(i/P_overload_sum) for i in P_overload] #give probabilities

        print("P_outaged is {},\n P_overload is {},\n".format(P_outaged,P_overload))


simulation(days, A, P_0, daily_lambda, max_F, F_0)
