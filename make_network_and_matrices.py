import numpy as np
from matplotlib import pyplot as plt


def make_network(n_layers):
    """
    Make a network with a certain number of layers. A layer is a cirkel around the previous layer with twice as many
    nodes.
    :param n_layers: Number of layers
    :return: range of nodes in each layer, list of lines
    """
    #  Compute the number of nodes per layer
    if n_layers == 1:
        nodes_per_layer = [1,]
    elif n_layers == 2:
        nodes_per_layer = [1, 3]
    elif n_layers > 2:
        nodes_per_layer = [1, 3]
        for i in range(2, n_layers):
            nodes_per_layer.append(nodes_per_layer[i-1]*2)
    else:
        nodes_per_layer = None
        exit('Error: n_layers should be an integer > 0')

    #  Compute which nodes are in each layer
    cum_nodes = np.cumsum(nodes_per_layer)
    cum_nodes = np.insert(cum_nodes, 0, 0)
    nodes_in_layer = []
    for i in range(n_layers):
        nodes_in_layer.append(list(cum_nodes[i:i+2]))

    #  Make connections between nodes
    #  First layer
    connections = [[0,1], [0,2], [0,3]]

    #  Additional layers
    for i in range(1, n_layers-1):
        # first and final node in current layer
        from_first = nodes_in_layer[i][0]
        from_final = nodes_in_layer[i][1]
        nodes_current = range(from_first, from_final)

        # first and final node in next layer
        to_first = nodes_in_layer[i+1][0]
        to_final = nodes_in_layer[i+1][1]
        nodes_next = range(to_first, to_final)

        for j in range(len(nodes_current)):
            connections.append([nodes_current[j], nodes_next[2*j]])
            connections.append([nodes_current[j], nodes_next[2*j + 1]])
    return nodes_in_layer, connections


def plot_network(buses_in_layer, lines, generator_indices=range(10, 22)):
    """
    Plots the network with lines in black, generators in red and loads in green
    :param buses_in_layer: range of nodes in each layer
    :param lines: list of lines
    :param generator_indices: id's of nodes which are a generator
    :return: None
    """
    x = np.array([])
    y = np.array([])
    for i in range(len(buses_in_layer)):
        b = range(buses_in_layer[i][0], buses_in_layer[i][1])
        n_buses = len(b)

        angles = np.linspace(0, 2.*np.pi, n_buses+1) + np.pi/n_buses
        angles = np.delete(angles, n_buses)
        x = np.append(x, np.cos(angles) * i)
        y = np.append(y, np.sin(angles) * i)

    for l in lines:
        x_line = [x[l[0]], x[l[1]]]
        y_line = [y[l[0]], y[l[1]]]
        plt.plot(x_line, y_line, c='black', zorder=0)

    plt.scatter(x, y, c='g', linewidths=5, zorder=1)
    plt.scatter(x[generator_indices], y[generator_indices], c='r', linewidths=5, zorder=2)
    plt.show()


def adjacency_matrix(connections):
    """
    Compute adjacency matrix
    :param connections: List of lines
    :return: Adjacency matrix M
    """
    n = np.amax(np.asarray(connections)) + 1
    M = np.zeros([n,n])
    for c in connections:
        M[c[0], c[1]] = 1
        M[c[1], c[0]] = 1
    return M


def matrix_a(connections, B):
    """
    Generates the matrix A given a list of lines and the corresponding susceptence matrix B
    :param connections: list of lines
    :param B: susceptence matrix
    :return: matrix A
    """
    m = len(connections)
    n = np.amax(np.asarray(connections)) + 1

    X = np.linalg.inv(B)
    N = np.zeros([m, n])

    for i, c in enumerate(connections):
        N[i, c[0]] = -B[c[0], c[1]]
        N[i, c[1]] = B[c[0], c[1]]

    A = N.dot(X)
    return A


buses_per_layer, lines = make_network(3)
M = adjacency_matrix(lines)

B = 3 * np.identity(M.shape[0]) - M
X = np.linalg.inv(B)
A = matrix_a(lines, B)

print(A)
plot_network(buses_per_layer, lines)
