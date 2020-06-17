import numpy as np

def power_vector(nr_buses, max_P_gen, NR_GENERATORS):

    # so that the max P for the generators is not exceeded (probably)
    high_P_load = (max_P_gen * (NR_GENERATORS - 1)) / (nr_buses) * 2

    buses_vect = []

    # TODO buses dict ws niet meer gebruiken
    buses = {}
    total_load = 0
    total_gen = 0

    # this is hardcoded for a system with 12 generators (depends on adjacency matrix)
    range_generators = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
    for i in range(1, nr_buses + 1):
        
        # this is hardcoded for a system with 12 generators (depends on adjacency matrix)
        if 10 < i < 23:
            buses[i] = {}
            buses[i]['type'] = 'generator'
        else:
            buses[i] = {}
            buses[i]['type'] = 'load'
            load = - np.random.uniform(high = high_P_load)
            buses_vect.append(load)
            total_load += load
                    
    for gen in range_generators:
        buses_vect = np.insert(buses_vect, gen, - (total_load / NR_GENERATORS))

    return buses_vect