import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import linprog
import random


class Network:
    def __init__(self, n_layers, generator_layer=0, generator_P_max=10, margin_F_max=2, labda=(1.00005, 0.25),
                 cost_weights=(1, 100), mu=1.005):
        """
        Make a network with a certain number of layers. A layer is a cirkel around the previous layer with twice as many
        nodes.

        PARAMETERS :

        n_layers: Number of layers
        generator_layer: Layer which contains the generators
        generator_P_max: Maximum power of one generator
        margin_F_max: F_max on day 0 is F_0 * margin_F_max
        labda: Daily power change, uniformly distributed between labda[0]/labda[1] and labda[0]*labda[1]
        cost_weights: Weight for generation and load shed
        mu: Line improvement factor

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
        self.P = self.set_p0(self.P_max)
        self.F = self.A.dot(self.P)
        self.F_max = np.abs(self.F * margin_F_max)
        self.M = np.abs(self.F / self.F_max)

        # Vector with cost function to redispatch power
        self.linprog_c = np.concatenate((cost_weights[0] * np.ones(len(self.generators)),
                                         -cost_weights[0] * np.ones(len(self.generators)),
                                         cost_weights[1] * np.ones(len(self.loads))),
                                        axis=0)

    ####################################################################################################################
    # Functions for initialization
    ####################################################################################################################
    def make_buses(self):
        '''
        Compute the number of nodes per layer

        RETURNS:

        n_nodes: Total number of nodes
        nodes_per_layer: List with the numbers of nodes per layer
        nodes_in_layer: List with the total number of nodes up to each layer

        '''

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
        '''
        Make lines between nodes

        RETURN:
        lines: List with tuples with the indices of the nodes that are connected with lines
        '''
        #First layer
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

        RETURN:
        M: adjacency matrix
        """
        M = np.zeros([self.n_nodes, self.n_nodes])
        for c in self.lines:
            M[c[0], c[1]] = 1
            M[c[1], c[0]] = 1
        return M

    def matrix_b(self, adjacency_matrix):
        """
        A matrix which its bij element is the susceptence of the transmission line joining bus i,j.

        RETURN:
        B: Matrix with susceptence of all lines = 1
        """
        diagonal = np.sum(adjacency_matrix, axis=1)
        B = np.diag(diagonal) - adjacency_matrix
        return B

    def b_with_outages(self, B, connected_nodes):
        """
        When lines break this affects B, this function returns B accounted for broken lines

        PARAMETERS:
        B: Previous matrix B
        connected_nodes: List containing a list with connected nodes for each line: [[0,1], [6,34],.. , [18, 28]]

        RETURN:
        new_B: Updated version of B where broken lines are taken into account
        """
        new_B = np.array(B)
        for nodes in connected_nodes:
            new_B[nodes[0], nodes[1]] = -0.0001
            new_B[nodes[1], nodes[0]] = -0.0001

        np.fill_diagonal(new_B, 0)
        diagonal = -np.sum(new_B, axis=1)
        new_B += np.diag(diagonal)
        return new_B

    def matrix_a(self, B):
        """
        A = NX
        N = the flow on the transmission line connecting bus i to j
        X = inv(B)

        PARAMETERS:
        B: susceptence matrix

        RETURN:
        A: matrix
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
        '''
        PARAMETERS:
        max_generator_power: maximum power limit

        RETURN:
        power: initial real power injection of a generator
        '''
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
        l = np.random.normal(self.labda[0], self.labda[1])
        self.P *= l
        self.P_max *= l

    def h0(self, overload_fraction):
        """
        Compute probability of initial failure
        PARAMETERS:
        overload_fraction:

        RETURN:

        """
        return 0.01 * overload_fraction**2

    def h1(self, overload_fraction):
        """
        Compute probability of failure due to overload
        """
        return 0.3 * overload_fraction**2

    def initial_failures(self):
        """
        Compute initial failures using h0
        RETURNS:
        line_indices: a list with the ids of failed lines
        connected_nodes: a list of coonected nodes of failed lines

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
        RETURNS:
        line_indices: list with the ids of failed lines
        connected_nodes: list with the connected nodes

        """
        line_indices = []
        connected_nodes = []
        for i, line in enumerate(self.lines):
            if self.M[i] > 0.95:
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

        PARAMETER:
        A: matrix

        RETURNS:
        sol.success: bool, true if the linprog algorith succeeds in finding an optimal solution
        load_shed: Flow trough lines

        """
        # Compute flows without redispatching power
        F = A.dot(self.P)

        # Equality constraint A_eq delta_P = b_eq: sum(delta_P) = 0
        A_eq = np.ones((1, len(self.linprog_c)))
        b_eq = np.array([0., ])

        # Inequality constraint: -Fmax - F < A delta_P < Fmax - F
        A_ub = np.concatenate((A[:, self.generators],
                               A[:, self.generators],
                               A[:, self.loads]),
                              axis=1)
        A_ub = np.concatenate((-A_ub, A_ub), axis=0)
        b_ub = np.concatenate((self.F_max + F, self.F_max - F))

        # Bounds for delta_P for generators (positive and negative) and loads
        bounds = []
        for generator in self.generators:
            bounds.append((0, self.P_max - self.P[generator]))
        for generator in self.generators:
            bounds.append((-self.P[generator], 0))
        for load in self.loads:
            bound = (0, -self.P[load])
            bounds.append(bound)

        # Solve new system
        sol = linprog(self.linprog_c, A_ub, b_ub, A_eq, b_eq, bounds=bounds, method='revised simplex')
        load_shed = 0

        if sol.success:
            n_generators = len(self.generators)
            delta_P = np.zeros(len(self.P))
            delta_P[self.generators] = sol.x[0:n_generators] + sol.x[n_generators:2 * n_generators]
            delta_P[self.loads] = sol.x[2 * n_generators::]
            load_shed += np.sum(sol.x[2 * n_generators::])
            print('Load shed:', load_shed)
            self.P = self.P + delta_P #hier moet 0 uitkomen als delta p == self.P

            # In this addition floating point errors occur if delta_P[i] = -self.P[i].
            # If this happens loads can have positive and generators negative power which causes errors.
            self.P[np.abs(self.P) < 1e-10] = 0

            self.F = A.dot(self.P)
            self.M = F / self.F_max
        return sol.success

    def improve_lines(self, failed_lines):
        """
        Improve lines which failed

        PARAMETER:
        failed_lines: the indices of failed lines

        """
        for line in failed_lines:
            self.F_max[line] *= self.mu

    def solar_panels(self):
        """
        Function to add solar energy to 20% of the loads and redistribute power of generators
        """
        max_P,sum_old_P = abs(max(self.P)),sum(self.P[self.loads])
        # print(random.sample(self.loads,round(0.2*len(self.loads))))
        solar_panels = [63, 31, 45, 76, 30, 72, 80, 50, 84, 54, 38, 19, 64, 87, 49, 77, 16, 82]
        for solar_panel in solar_panels:
            self.load_efficiency = random.random()
            self.P[solar_panel] = (max_P * self.load_efficiency)+self.P[solar_panel]
        sum_new_P = sum(self.P[self.loads])
        solar_energy = abs(sum_new_P-sum_old_P)
        subtract_P = solar_energy/len(self.generators)
        for generator in self.generators:
            self.P[generator] -= subtract_P
    ####################################################################################################################
    # Function to simulate a day
    ####################################################################################################################
    def simulate_day(self):
        """
        Function to simulate a day

        RETURNS:
        failed_lines: The ids of failed lines at the end of the day
        load_shed_today: The shedding load of the day
        """

        # Copy of adjacency matrix where outages lines will be removed

        # Slow dynamics
        self.update_P()
        self.F = self.A.dot(self.P)
        self.M = self.F / self.F_max

        failed_lines = []
        line_ids, linked_nodes = self.initial_failures()
        B = np.array(self.B)
        while len(line_ids) > 0:
            failed_lines.append(line_ids)

            B = self.b_with_outages(B, linked_nodes)
            A = self.matrix_a(B)

            success = self.redispatch_power(A)

            if success:
                line_ids, linked_nodes = self.overload_failures()
            else:
                line_ids = []

        self.improve_lines(failed_lines)
        return failed_lines

    ####################################################################################################################
    # Functions for visualization and representation
    ####################################################################################################################
    def report_params(self):
        '''
        RETURNS:
        The number of nodes
        The number of generators
        THe number of loads
        The power levels
        The flow
        The overload fraction
        '''
        print('Number of nodes:', self.n_nodes)
        print('Number of generators:', len(self.generators))
        print('Number of loads:', len(self.loads))
        print('Power levels:', self.P)
        print('Flow:', self.F)
        print('Overload fraction:', self.M)

    def plot_network(self):
        """
        Plots the network with lines in black, generators in red and loads in green

        RETURN:
        Figure with the network
        """

        x = np.array([])
        y = np.array([])
        for i in range(len(self.nodes_in_layer)):
            b = range(self.nodes_in_layer[i][0], self.nodes_in_layer[i][1])
            n_buses = len(b)

            angles = np.linspace(0, 2.*np.pi, n_buses+1) + np.pi/n_buses
            angles = np.delete(angles, n_buses)
            x = np.append(x, np.cos(angles) * i)
            y = np.append(y, np.sin(angles) * i)

        for l in self.lines:
            x_line = [x[l[0]], x[l[1]]]
            y_line = [y[l[0]], y[l[1]]]
            plt.plot(x_line, y_line, c='black', zorder=0)

        plt.scatter(x, y, c='g', linewidths=5, zorder=1)
        plt.scatter(x[self.generators], y[self.generators], c='r', linewidths=5, zorder=2)
        plt.show()

N = Network(6, 2)

N.solar_panels()
for i in range(100):
    print('day', i)
    lines = N.simulate_day()
    print('failed lines:', lines)
B = N.b_with_outages(N.B, [])
