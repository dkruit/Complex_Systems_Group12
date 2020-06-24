import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import linprog
from scipy.interpolate import UnivariateSpline
import random
import os
import csv

N_LAYERS = 6
GEN_LAYER = 2
N_SIMS = 1
N_DAYS = 20
NETWORK_FIGS = False
PLOT_FIGS = True
EXPERIMENT_NAME = "test_14u38"

class Network:
    def __init__(self, n_layers, generator_layer=0, generator_P_max=10, margin_F_max=1.2, labda=(1.00005, 1.005),
                 cost_weights=(1, 100), mu=1.005):
        """
        Make a network with a certain number of layers. A layer is a cirkel around the previous layer with twice as many
        nodes.
        :param n_layers: Number of layers
        :param generator_layer: Layer which contains the generators
        :param generator_P_max: Maximum power of one generator
        :param margin_F_max: F_max on day 0 is F_0 * margin_F_max
        :param labda: Daily power change is a ubiform distribution between labda[0]/labda[1] and labda[0]*labda[1]
        :param cost_weights: Weight for generation and load shed
        """
        # Set parameters
        self.labda = labda
        self.mu = mu

        # Initialize network
        self.n_layers = n_layers
        self.n_nodes, self.nodes_per_layer, self.nodes_in_layer = self.make_buses()
        self.generators = range(self.nodes_in_layer[generator_layer][0], self.nodes_in_layer[generator_layer][1])
        self.loads = list(range(self.generators[0])) + list(range(self.generators[-1]+1, self.n_nodes))
        self.lines = self.make_lines()

        # Initialize matrices
        self.M_adj = self.adjacency_matrix()
        self.B = self.matrix_b(self.M_adj)
        self.A = self.matrix_a(self.B)

        # Initialize P, F, Fmax and M
        self.P_max = generator_P_max
        self.P_slow = self.set_p0(self.P_max)
        self.P_fast = np.array(self.P_slow)

        self.F = self.A.dot(self.P_slow)
        self.F_max = np.abs(self.F * margin_F_max)
        self.F_max_fast = np.array(self.F_max)
        self.M = np.abs(self.F / self.F_max)

        # Vector with cost function to redispatch power
        self.linprog_c = np.concatenate((cost_weights[0] * np.ones(len(self.generators)),
                                         -cost_weights[0] * np.ones(len(self.generators)),
                                         cost_weights[1] * np.ones(len(self.loads))),
                                        axis=0)

        # Vector to store results
        self.load_shed = []
        self.n_failed = []

        # Storage for analysis and visualisation
        self.failures_store = []
        self.maximal_F = np.max(self.F_max)
        self.day = 0
        self.store_F = []
        self.store_P = []
        self.store_P_redispatched = []
        self.delta_P = np.zeros(len(self.P_fast))

        # Vector to store results
        self.load_shed = []
        self.n_failed = []

    ####################################################################################################################
    # Functions for initialization
    ####################################################################################################################
    def make_buses(self):
        #  Compute the number of nodes per layer
        if self.n_layers == 1:
            nodes_per_layer = [1, ]
        elif self.n_layers == 2:
            nodes_per_layer = [1, 3]
        elif self.n_layers > 2:
            nodes_per_layer = [1, 3]
            for i in range(2, self.n_layers):
                nodes_per_layer.append(nodes_per_layer[i - 1] * 2)
        else:
            nodes_per_layer = None
            exit('Error: self.n_layers should be an integer > 0')

        #  Compute which nodes are in each layer
        n_nodes = np.sum(nodes_per_layer)
        cum_nodes = np.cumsum(nodes_per_layer)
        cum_nodes = np.insert(cum_nodes, 0, 0)
        nodes_in_layer = []
        for i in range(self.n_layers):
            nodes_in_layer.append(list(cum_nodes[i:i + 2]))
        return n_nodes, nodes_per_layer, nodes_in_layer

    def make_lines(self):
        #  Make lines between nodes
        #  First layer
        lines = [[0, 1], [0, 2], [0, 3]]

        #  Additional layers
        for i in range(1, self.n_layers - 1):
            # first and final node in current layer
            from_first = self.nodes_in_layer[i][0]
            from_final = self.nodes_in_layer[i][1]
            nodes_current = range(from_first, from_final)

            # first and final node in next layer
            to_first = self.nodes_in_layer[i + 1][0]
            to_final = self.nodes_in_layer[i + 1][1]
            nodes_next = range(to_first, to_final)

            for j in range(len(nodes_current)):
                lines.append([nodes_current[j], nodes_next[2 * j]])
                lines.append([nodes_current[j], nodes_next[2 * j + 1]])
        return lines

    def adjacency_matrix(self):
        """
        Compute adjacency matrix
        """
        M = np.zeros([self.n_nodes, self.n_nodes])
        for c in self.lines:
            M[c[0], c[1]] = 1
            M[c[1], c[0]] = 1
        return M

    def matrix_b(self, adjacency_matrix):
        """
        Computes the adjacency matrix. For now uses susceptence of 1 for all lines.
        :return:
        """
        diagonal = np.sum(adjacency_matrix, axis=1)
        B = np.diag(diagonal) - adjacency_matrix
        return B

    def b_with_outages(self, B, connected_nodes):
        """
        When lines break this affects B, this function returns B accounted for broken lines
        :param B: Previous matrix B
        :param connected_nodes: List containing a list with connected nodes for each line: [[0,1], [6,34],.. , [18, 28]]
        :return: Updated version of B where broken lines are taken into account
        """
        new_B = np.array(B)
        for nodes in connected_nodes:
            new_B[nodes[0], nodes[1]] = -0.0000001
            new_B[nodes[1], nodes[0]] = -0.0000001

        np.fill_diagonal(new_B, 0)
        diagonal = -np.sum(new_B, axis=1)
        new_B += np.diag(diagonal)
        return new_B

    def matrix_a(self, B):
        """
        Generates the matrix A given a list of lines and the corresponding susceptence matrix B
        :return: matrix A
        """
        m = len(self.lines)

        X = np.linalg.inv(B)
        N = np.zeros([m, self.n_nodes])

        for i, c in enumerate(self.lines):
            N[i, c[0]] = -B[c[0], c[1]]
            N[i, c[1]] = B[c[0], c[1]]

        A = N.dot(X)
        return A

    def set_p0(self, max_generator_power):
        n_generators = len(self.generators)
        n_loads = self.n_nodes - n_generators
        max_load_power = max_generator_power * n_generators / n_loads

        power = np.zeros(self.n_nodes)
        total_load = 0

        for l in self.loads:
            power[l] = -np.random.uniform(high=max_load_power)
            total_load += power[l]

        generator_power = -total_load / n_generators
        power[self.generators] = generator_power
        return power

    ####################################################################################################################
    # Functions required for simulation
    ####################################################################################################################
    def update_P(self):
        """
        Update P using lambda
        """
        l = np.random.uniform(self.labda[0], self.labda[0]*self.labda[1]) ** np.random.choice([-1, 1])
        print('l', l)
        self.P_slow *= l
        self.P_max *= l

    def h0(self, overload_fraction):
        """
        Compute probability of initial failure
        """
        return 0.01 * overload_fraction**2

    def h1(self, overload_fraction):
        """
        Compute probability of failure due to overload
        """
        return 0.1 * overload_fraction**2

    def initial_failures(self):
        """
        Compute initial failures using h0
        :return: List of ids and connected nodes of failed lines
        """
        line_indices = []
        connected_nodes = []
        for i, line in enumerate(self.lines):
            p = self.h0(self.M[i])
            a = np.random.uniform(0, 1)
            if p > a:
                line_indices.append(i)
                connected_nodes.append(line)
        return line_indices, connected_nodes

    def overload_failures(self):
        """
        Compute secondary failures due to overload
        :return: List of failed lines
        """
        line_indices = []
        connected_nodes = []
        for i, line in enumerate(self.lines):
            if self.M[i] > 0.99 and self.F_max_fast[i] > 0.000001:
                p = self.h1(self.M[i])
                if p > np.random.uniform(0, 1):
                    line_indices.append(i)
                    connected_nodes.append(line)
        return line_indices, connected_nodes
    pass

    def redispatch_power(self, A):
        """
        Find new solution to redispatch power after failures occur
        Delta_P for the generators is placed in front of loads split in positive and negative part:
        delta_P = [gen_1+, gen_2+ ... gen_n+, gen_1-, gen_2- ... gen_n-, load_1, load_2 ... load_m] for n generators
        and m loads.
        :return: Flow trough lines
        """
        # Compute flows without redispatching power
        F = A.dot(self.P_fast)

        # Equality constraint A_eq delta_P = b_eq: sum(delta_P) = 0
        A_eq = np.ones((1, len(self.linprog_c)))
        b_eq = np.array([0., ])

        # Inequality constraint: -Fmax - F < A delta_P < Fmax - F
        A_ub = np.concatenate((A[:, self.generators],
                               A[:, self.generators],
                               A[:, self.loads]),
                              axis=1)
        A_ub = np.concatenate((-A_ub, A_ub), axis=0)
        b_ub = np.concatenate((self.F_max_fast + F, self.F_max_fast - F))

        # Bounds for delta_P for generators (positive and negative) and loads
        bounds = []
        for generator in self.generators:
            bound = self.P_max - self.P_fast[generator]
            # Due to floating point error this may become negative.
            if bound < 0:
                bound = 0.
            bounds.append((0, bound))
        for generator in self.generators:
            bounds.append((-self.P_fast[generator], 0))
        for load in self.loads:
            bound = (0, -self.P_fast[load])
            bounds.append(bound)

        # Solve new system
        sol = linprog(self.linprog_c, A_ub, b_ub, A_eq, b_eq, bounds=bounds, method='revised simplex')
        if sol.success:
            n_generators = len(self.generators)
            self.delta_P = np.zeros(len(self.P_fast))
            self.delta_P[self.generators] = sol.x[0:n_generators] + sol.x[n_generators:2 * n_generators]
            self.delta_P[self.loads] = sol.x[2 * n_generators::]
            load_shed = np.sum(self.delta_P[self.loads])
            self.P_fast += self.delta_P

            # In this addition floating point errors occur if delta_P[i] = -self.P[i].
            # If this happens loads can have positive and generators negative power which causes errors.
            self.P_fast[np.abs(self.P_fast) < 1e-10] = 0

            self.F = A.dot(self.P_fast)
            self.M = np.abs(self.F) / self.F_max
        else:
            load_shed = 0
        return sol.success, load_shed

    def improve_lines(self, failed_lines):
        """
        Improve lines which failed
        """
        for line in failed_lines:
            self.F_max[line] *= self.mu

    def solar_panels(self):
        """
        Function to add solar energy to 20% of the loads
        """

        max_P,sum_old_P = abs(max(self.P_slow)),sum(self.P_slow[self.loads])
        solar_panels = random.sample(self.loads,round(0.2*len(self.loads)))
        for i in solar_panels:
            self.load_efficiency = random.random()
            self.P_slow[i] = (max_P * self.load_efficiency)+self.P_slow[i]
        sum_new_P = sum(self.P_slow[self.loads])
        solar_energy = abs(sum_new_P-sum_old_P)
        subtract_P = solar_energy/len(self.generators)
        for i in self.generators:
            self.P_slow[i] -= subtract_P

    ####################################################################################################################
    # Function to simulate a day
    ####################################################################################################################
    def simulate_day(self):
        """
        Function to simulate a day
        """
        # Copy of adjacency matrix where outages lines will be removed

        # Slow dynamics
        self.update_P()
        self.F = self.A.dot(self.P_slow)
        self.F_max_fast = np.array(self.F_max)
        self.M = np.abs(self.F) / self.F_max
        self.P_fast = np.array(self.P_slow)
        self.plot_network(function = "process", title = "a_updated_p")
        # print(np.mean(np.abs(self.P_slow)), np.mean(self.M))

        failed_lines = []
        load_shed_today = 0.
        line_ids, linked_nodes = self.initial_failures()
        B = np.array(self.B)
        self.plot_network(function = "process", title = "b_initial_failure", failed_lines_plot = line_ids)

        while len(line_ids) > 0:
            for line_id in line_ids:
                if line_id not in failed_lines:
                    failed_lines.append(line_id)
            self.F_max_fast[line_ids] = 0.

            """
            B = self.b_with_outages(B, linked_nodes)
            try:
                A = self.matrix_a(B)
            except:
                print('Singular matrix, ending day')
                success = False
            else:
                success, load_shed = self.redispatch_power(A)
                load_shed_today += load_shed
            if success:
                line_ids, linked_nodes = self.overload_failures()
            else:
                line_ids = []
            """
            success, load_shed = self.redispatch_power(self.A)
            load_shed_today += load_shed

            if success:
                line_ids, linked_nodes = self.overload_failures()
            else:
                line_ids = []

        self.improve_lines(failed_lines)
        if len(failed_lines) > 0:
            self.load_shed.append(load_shed_today / np.sum(self.P_slow[self.generators]))
            self.n_failed.append(len(failed_lines))
        self.plot_network(function = "process", title = "c_after_redispatch", failed_lines_plot = failed_lines)
        self.day += 1
        return failed_lines, load_shed_today

    ####################################################################################################################
    # Functions for visualization and representation
    ####################################################################################################################
    def report_params(self):
        print('Number of nodes:', self.n_nodes)
        print('Number of generators:', len(self.generators))
        print('Number of loads:', len(self.loads))
        print('Power levels:', self.P_slow)
        print('Flow:', self.F)
        print('Overload fraction:', self.M)

    def plot_network(self, function = "network", title = "", failed_lines_plot = [], failed_nodes_plot = []):
        """
        Plots the network with lines in black, generators in red and loads in green
        """
        # x = np.array([])
        # y = np.array([])
        # for i in range(len(self.nodes_in_layer)):
        #     b = range(self.nodes_in_layer[i][0], self.nodes_in_layer[i][1])
        #     n_buses = len(b)

        #     angles = np.linspace(0, 2.*np.pi, n_buses+1) + np.pi/n_buses
        #     angles = np.delete(angles, n_buses)
        #     x = np.append(x, np.cos(angles) * i)
        #     y = np.append(y, np.sin(angles) * i)

        # for l in self.lines:
        #     x_line = [x[l[0]], x[l[1]]]
        #     y_line = [y[l[0]], y[l[1]]]
        #     plt.plot(x_line, y_line, c='black', zorder=0)

        # plt.scatter(x, y, c='g', linewidths=5, zorder=1)
        # plt.scatter(x[self.generators], y[self.generators], c='r', linewidths=5, zorder=2)
        # plt.show()
        if NETWORK_FIGS:
            x = np.array([])
            y = np.array([])
            for i in range(len(self.nodes_in_layer)):
                b = range(self.nodes_in_layer[i][0], self.nodes_in_layer[i][1])
                n_buses = len(b)

                angles = np.linspace(0, 2.*np.pi, n_buses+1) + np.pi/n_buses
                angles = np.delete(angles, n_buses)
                x = np.append(x, np.cos(angles) * i)
                y = np.append(y, np.sin(angles) * i)

                # # todo
                # if self.delta_P[i] > 0 and self.P_fast[i] < 0 and self.P_fast[i] + self.delta_P[i] < 0.005:

                #     # this should be black, because
                #     plt.scatter(x, y, c='k', edgecolors = 'k', linewidths=1, zorder=1, s = 120)
                # else:
                #     plt.scatter(x, y, c='g', edgecolors = 'g', linewidths=1, zorder=1, s = 120)

            for i, l in enumerate(self.lines):
                x_line = [x[l[0]], x[l[1]]]
                y_line = [y[l[0]], y[l[1]]]
                if function == "network":
                    plt.plot(x_line, y_line, c='black', zorder=0)
                else:
                    if i in failed_lines_plot:
                        plt.plot(x_line, y_line, 'k--', zorder=0, linewidth = 1 + 2 * (self.F_max[i] / self.maximal_F))
                    else:
                        plt.plot(x_line, y_line, c = self.get_linecolor(i), zorder=0, linewidth = 1 + 2 * (self.F_max[i] / self.maximal_F))

            # hier die forloop zetten
            # plt.scatter(x, y, c='g', edgecolors = 'g', linewidths=1, zorder=1, s = 120)
            plt.scatter(x[self.generators], y[self.generators], c='w', edgecolors = 'g', linewidths=1, zorder=2, s = 120)
            plt.title(label = title + " day " + str(self.day))
            plt.savefig(EXPERIMENT_NAME + "/day" + str(self.day) + title + ".jpg")
            plt.close()

    def get_linecolor(self, line_number):
        # colors0 = cm.get_cmap("Greens")
        # colors1 = cm.get_cmap("Reds")
        overflow = self.M[line_number]
        if overflow <= 0.99:
            return "g"
        else:
            return "r"





# N = Network(6, 2)
# # N.report_params()
# # N.solar_panels()
# for i in range(5000):
#     print('\nday', i)
#     lines, load_shed = N.simulate_day()
#     print('failed lines:', lines)
#     print('load shed:', load_shed)
#     print('power', -np.sum(N.P_slow[N.loads]))
#     print(np.mean(np.abs(N.P_slow)), np.mean(N.M))
#     # if a == 1:
#         # break

# plt.hist(N.n_failed)
# plt.show()


def start_simulation(nr_layers, generator_layer, n_simulations, n_days):
    if not os.path.exists(EXPERIMENT_NAME):
        os.makedirs(EXPERIMENT_NAME)
    for i in range(n_simulations):
        all_failed_lines = []
        shedded_loads = []
        print('\nsimulation', i)
        N = Network(nr_layers, generator_layer)
        for j in range(n_days):
            lines, load_shed = N.simulate_day()
            if len(lines) > 0:
                all_failed_lines.append(len(lines))
            if load_shed > 0:
                shedded_loads.append(load_shed / np.sum(N.P_slow[N.loads]))
            shedded_loads.append(load_shed)
            print('failed lines:', lines)
            print('load shed:', load_shed)
            print(np.mean(np.abs(N.P_slow)), np.mean(N.M))

        B = N.b_with_outages(N.B, [])
        failed_lines_file = open(EXPERIMENT_NAME + "/flines_simulation_" + str(i) + ".txt", 'w')
        failed_lines_file.write(str(all_failed_lines))
        failed_lines_file.close()
        shedded_loads_file = open(EXPERIMENT_NAME + "/shedded_simulation_" + str(i) + ".txt", 'w')
        shedded_loads_file.write(str(shedded_loads))
        shedded_loads_file.close()

    if PLOT_FIGS:

        # hist, bins, _ = plt.hist(all_failed_lines)

        # # histogram on log scale.
        # # Use non-equal bin sizes, such that they look equal on log scale.
        # logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
        # plt.show()
        plt.hist(all_failed_lines, log = True)
        plt.xlabel(xlabel = "Nr of failed lines")
        plt.ylabel(ylabel = "Frequency")
        plt.savefig(EXPERIMENT_NAME + "/frequency_failedlines_hist.jpg")
        plt.show()

        line_outages = [0] * 94
        for nr in all_failed_lines:
            line_outages[nr] += 1

        plt.plot(range(0, len(line_outages)), line_outages)
        plt.xlabel(xlabel = 'Nr of failed lines')
        plt.ylabel(ylabel = 'Nr of events')
        plt.grid(b = True)
        plt.savefig(EXPERIMENT_NAME + "/Nrevents_failedlines.jpg")
        plt.show()

        x_outages = []
        y_outages = []

        for i, value in enumerate(line_outages):
            if value > 0:
                x_outages.append(i)
                y_outages.append(value)

        for i, value in enumerate(y_outages):
            y_outages[i] = value/len(all_failed_lines)

        plt.plot(x_outages, y_outages)
        plt.xlabel(xlabel = "Nr of failed lines")
        plt.ylabel(ylabel = "Probability")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b = True)
        plt.savefig(EXPERIMENT_NAME + "/probability_failedlines_log.jpg")
        plt.show()

        plt.hist(shedded_loads, density = True)
        plt.show()

        n = 20
        p, x = np.histogram(N.load_shed, bins=np.logspace(-2, 0, n))
        x = x[:-1] + (x[1] - x[0])/2   # convert bin edges to centers
        f = UnivariateSpline(x, p, s=n)
        plt.plot(x, f(x))
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(xlabel = "Load shed")
        plt.ylabel(ylabel = "Frequency")
        plt.xticks((1e-2, 1e-1, 1))
        plt.grid(b = True)
        plt.savefig(EXPERIMENT_NAME + "/frequency_loadshed_log.jpg")
        plt.show()



start_simulation(N_LAYERS, GEN_LAYER, N_SIMS, N_DAYS)