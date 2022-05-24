import networkx as nx
import numpy as np
import random
from tqdm import tqdm


def read_network_static_singular(file, comment="%", delim="\t"):
    g = nx.MultiDiGraph() # directed graph with multiple edges

    with open(file, "r") as f:
        # iterate over the file lines
        for l in tqdm(f):
            # Skip comment lines
            if l[0] != comment:

                # put edge data into a list of integers
                data = list(map(int, l.strip().split(delim)))

                # always add the edge between the nodes, even if already existing
                g.add_edge(data[0], data[1])

    # Compute the weight based on the number of nodes
    weights = dict()

    for e in g.edges():
        if not e in weights.keys():
            weights[e] = 1
        else:
            weights[e] += 1

    # Finally, creating the Non-Multi Directed Graph (no multiple edges)
    gf = nx.DiGraph()

    for k in weights.keys():
        gf.add_edge(k[0], k[1], weight=weights[k])

    return gf



def get_biggest_weakly_connected_component(g):
    biggest_comp_nodes = max(nx.weakly_connected_components(g), key=len)
    return g.subgraph(biggest_comp_nodes)



def read_network_temporal(file, comment="%", delim="\t"):
    time_dict = {} # initialize dictionary
    nodes = set() # set of network nodes

    # open the file in read mode
    with open(file, "r") as f:

        # iterate over the lines
        for l in f:
            # discard comment lines
            if l[0] != comment:

                # put edge data into a list of integers
                data = list(map(int, l.strip().split(delim)))

                # add nodes to the set of nodes
                nodes.update({data[0], data[1]})

                # new timestamp
                if not data[3] in time_dict.keys():

                    # add edge as a tuple in the list of the corresponding timestamp
                    time_dict[data[3]] = [(data[0], data[1])]

                # already present timestamp/key
                else:
                    # add new tuple/edge to the list indexed by that timestamp
                    time_dict[data[3]] += [(data[0], data[1])]

    # returns the dictionary and the nodes set
    return time_dict, nodes



# immunized_nodes = list of (eventually) immune nodes
# infection_prob = p parameter of infecting a susceptible node
# time_dict = dictionary sorted by timestamps
# seed = initial infected node/nodes

def simulate_SI(time_dict, seed, nodes, n_nodes, infection_prob, immunized_nodes= []):
    """
    immunized_nodes = list of (eventually) immune nodes
    infection_prob = p parameter of infecting a susceptible node
    time_dict = dictionary sorted by timestamps
    seed = initial infected node/nodes
    """
    # initialize infection_time of each node to 'inf'
    infection_times = {} # dictionary with infection time of each node
    for node in sorted(nodes):
        #infection_times[node] = float('inf')
        infection_times[node] = np.inf

    # initialize the infection time of the seed nodes to the first timestamp
    for el in seed: # N.B seed needs to be a list
        infection_times[el] = list(time_dict.keys())[0]

    # infected nodes count
    infected_nodes_count = len(seed)

    # rho (percentage of infected nodes) at each timestamp
    # initialize with zeros
    rho_list = np.zeros(len(time_dict)) # Initial rho + one rho per timestep

    # initial rho
    #rho_list[0] = infected_nodes_count

    # iterate over dictionary following ordered timestamps
    for i, t in enumerate(time_dict.keys()): # notice that t is the current timestamp

        # loop through all the edges
        for edge in time_dict[t]:

            # check if first node in the edge is infected at current timestamp
            if infection_times[edge[0]] <= t:

                # destination node not already infected
                if infection_times[edge[1]] > t:

                    # generate random number in [0,1]
                    prob = random.random()

                    # infect destination node and set its time of infection
                    if prob < infection_prob:
                        infection_times[edge[1]] = t

                        # increase infected nodes count
                        infected_nodes_count += 1

        # after each timestamp, update rho list
        rho_list[i] = infected_nodes_count

    # return tho and nodes infection time
    return rho_list/n_nodes, infection_times



def run_multiple_simulations(n_simulations, time_dict, seed, nodes, n_nodes, infection_prob, immunized_nodes= []):
    rho_lists = []
    #rho_list = np.zeros(len(time_dict)) # Initial rho + one rho per timestep
    infection_times_lists = []

    # main loop
    for i in tqdm(range(n_simulations)):
        rho_list, infection_times = simulate_SI(time_dict, seed, nodes, n_nodes, infection_prob, immunized_nodes)
        rho_lists.append(rho_list)
        infection_times_lists.append(infection_times)

    # we need to return now the average rho and infection time over the lists
    rho_final = np.zeros(len(time_dict))
    for i in range(n_simulations):
        for j in range(len(rho_final)):
            rho_final[j] +=  rho_lists[i][j]

    return rho_final/n_simulations