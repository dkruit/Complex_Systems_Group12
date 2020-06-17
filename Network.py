import numpy as np
from matplotlib import pyplot as plt


class Network:
    def __init__(self, n_layers, generator_layer=0, generator_P_max=10, margin_F_max=1.05, labda=(1.0005, 0.0005)):
        """
        Make a network with a certain number of layers. A layer is a cirkel around the previous layer with twice as many
        nodes.
        :param n_layers: Number of layers
        :param generator_layer: Layer which contains the generators
        :param generator_P_max: Maximum power of one generator
        :param margin_F_max: F_max on day 0 is F_0 * margin_F_max
        :param labda: Mean and standard deviation of average daily power change
        """
        # Set parameters
        self.labda = labda

        # Initialize network
        self.n_layers = n_layers
        self.n_nodes, self.nodes_per_layer, self.nodes_in_layer = self.make_buses()
        self.generators = range(self.nodes_in_layer[generator_layer][0], self.nodes_in_layer[generator_layer][1])
        self.loads = list(range(self.generators[0])) + list(range(self.generators[-1]+1, self.n_nodes))
        self.lines = self.make_lines()

        # Initialize matrices
        self.M_adj = self.adjacency_matrix()
        self.B = self.matrix_b()
        self.A = self.matrix_a()

        # Initialize P, F, Fmax and M
        self.P_max = generator_P_max
        self.P = self.set_p0(self.P_max)
        self.F = self.A.dot(self.P)
        self.F_max = np.abs(self.F * margin_F_max)
        self.M = np.abs(self.F / self.F_max)

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

    def matrix_b(self):
        """
        Computes the adjacency matrix. For now uses susceptence of 1 for all lines.
        :return:
        """
        diagonal = np.sum(self.M_adj, axis=1)
        B = np.diag(diagonal) - self.M_adj
        return B

    def matrix_a(self):
        """
        Generates the matrix A given a list of lines and the corresponding susceptence matrix B
        :return: matrix A
        """
        m = len(self.lines)

        X = np.linalg.inv(self.B)
        N = np.zeros([m, self.n_nodes])

        for i, c in enumerate(self.lines):
            N[i, c[0]] = -self.B[c[0], c[1]]
            N[i, c[1]] = self.B[c[0], c[1]]

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
        l = np.random.normal(self.labda[0], self.labda[1])
        self.P *= l
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
        return 0.5 * overload_fraction**2

    def initial_failures(self):
        """
        Compute initial failures using h0
        :return: List of ids and connected nodes of failed lines
        """
        line_indices = []
        connected_nodes = []
        for i, line in enumerate(self.lines):
            p = self.h0(self.M[i])
            if p > np.random.uniform(0, 1):
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
            p = self.h1(self.M[i])
            if p > np.random.uniform(0, 1):
                line_indices.append(i)
                connected_nodes.append(line)
        return line_indices, connected_nodes
    pass

    def redispatch_power(self):
        """
        Find new solution to redispatch power after failures occur
        :return: List of failed lines
        """
        pass

    def improve_lines(self):
        """
        Improve lines which failed
        """
        pass

    ####################################################################################################################
    # Function to simulate a day
    ####################################################################################################################
    def simulate_day(self):
        """
        Function to simulate a day
        """

        # Slow dynamics
        self.update_P()
        self.F = self.A.dot(self.P)
        self.M = np.max(self.F / self.F_max)

        self.initial_failures()
        self.redispatch_power()

        self.improve_lines()

    ####################################################################################################################
    # Functions for visualization
    ####################################################################################################################
    def plot_network(self):
        """
        Plots the network with lines in black, generators in red and loads in green
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


N = Network(6, 3)
for i in range(10):
    a, b = N.initial_failures()
    print(a, b)
N.plot_network()
