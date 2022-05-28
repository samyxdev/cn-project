import networkx as nx
import numpy as np
import random
from tqdm import tqdm
import matplotlib.pyplot as plt

# for plots
plt.rcParams["figure.figsize"] = (12, 6)
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams.update({'font.size': 20})
plt.style.use('seaborn-darkgrid')


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

def get_biggest_strongly_connected_component(g):
    biggest_comp_nodes = max(nx.strongly_connected_components(g), key=len)
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

    # returns the dictionary(sorted by timestamps) and the order nodes list
    return dict(sorted(td.items())), sorted(nodes)

def read_network_temporal_filtered(file, nodes_to_keep, comment="%", delim="\t"):
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

                # Because we have the biggest weakly compnent graph, we just
                # need to to check one side of the edge to ensure that this is a
                # node to skip (or to keep)
                if not data[0] not in nodes_to_keep:
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

    # returns the dictionary(sorted by timestamps) and the order nodes list
    return dict(sorted(time_dict.items())), sorted(nodes)



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
    for node in nodes:
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

    # build average rho and infection times over the lists
    rho_final = np.zeros(len(time_dict))

    infection_times = {node_id:0 for node_id in nodes}
    for i in range(n_simulations):
        for j in range(len(rho_final)):
            rho_final[j] +=  rho_lists[i][j]

        #for j in range(n_nodes):
        #    infection_times[i]


    return rho_final/n_simulations


def run_multiple_simulations2(n_simulations, time_dict, seed, nodes, n_nodes, infection_prob, immunized_nodes= []):
    """Apply multiple simulations having defined seed nodes, infection probability and immunized_nodes if any.

    Return the average infection coefficient (rho) and average infection time, over all simulations"""
    rho_lists = []
    #rho_list = np.zeros(len(time_dict)) # Initial rho + one rho per timestep
    infection_times_lists = []


    # main loop
    for i in tqdm(range(n_simulations)):
        rho_list, infection_times = simulate_SI(time_dict, seed, nodes, n_nodes, infection_prob, immunized_nodes)
        rho_lists.append(rho_list)
        infection_times_lists.append(infection_times.values())

    # build average rho and infection times over the lists
    rho_final = np.zeros(len(time_dict))
    infection_times_final = np.zeros(len(nodes))
    for i in range(n_simulations):
        for j in range(len(rho_final)):
            rho_final[j] +=  rho_lists[i][j]
        for j in range(n_nodes):
            infection_times_final[j] += infection_times_lists[i][j]


    return rho_final/n_simulations, infection_times_final/n_simulations


# prevalence plot
def plot_prevalence(time_list, rho_lists, labels=[]):
    plt.figure(figsize=(12,6))
    for i in range(len(rho_lists)):
        plt.plot(time_list, rho_lists[i], label = labels[i])
        plt.title("ρ in function of time", fontsize=15)
        plt.xlabel("Time after first event (from zero)", fontsize=15)
        plt.ylabel("ρ", fontsize=15)
        plt.yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        plt.legend(loc='upper left')
        plt.style.use('seaborn-darkgrid')
        plt.tight_layout()


def random_group(possible_states, probabilites):
        """Randomly return a certain state from possible_states
        based on a set of probabilities of which the sum must be 1
        """
        assert sum(probabilites) == 1.0

        choices = []

        for i, st in enumerate(possible_states):
            choices += [st]*int(probabilites[i] * 10)

        return np.random.choice(choices)





def simulation_step(tx_state, rx_state, rx_inf_prob):
        """Apply one step of the simulation to the nodes, based on their state (Infected or not)
        and on the reciever infection probability rx_inf.
        Return the new state and a int equals to 1 if the state has changed, 0 otherwise

        TODO: In case of implemting SI with Recovery (SIR), then choice must be -1 if recovery
        and 1 for infection (for simulation_loop)
        """
        if tx_state and not rx_state:
            choice = np.random.binomial(1, rx_inf_prob)
            return bool(choice), choice
        else:
            return rx_state, False

        #return np.random.binomial(1, rx_inf_prob) if tx_state and not rx_state else rx_state


def simulation_loop(time_dict, states, n_nodes, inf_prob):
        """Applies an SI simulation based on the contact made at each timestamp.
        Returns an array containing the evolution of the rho factor at each timestamp.

        TODO: Add the support for different probabilities depending on the group of the node
        """
        rho_list = np.zeros(1+len(time_dict)) # Initial rho + one rho per timestep

        for i, t in enumerate(time_dict.keys()):
            # For each contact

            # To speed up rho calculations, we'll keep track of the number of infected by monitoring
            # if the state of nodes changed. If they went to infected, we add one to the previously
            # saved rho.
            temp_rho = 0
            for contact in time_dict[t]:
                states[contact[1]], changed = simulation_step(states[contact[0]], states[contact[1]], inf_prob)

                temp_rho += changed

            # We always sums the changed as there is no recovery (cf. TODO about SIR in simulation_step)
            rho_list[i+1] = rho_list[i] + temp_rho

        return rho_list/(n_nodes+1)

'''
# Simple case: only one group, random selection of initial infected (0 or 1)
states = dict(zip(nodes, np.random.binomial(1, 0.05, n_nodes)))

# More complicated : Case with different groups
grps = (1, 2, 3)
probs = (0.1, 0.4, 0.5)
#states = dict(zip(nodes, [random_group(grps, probs) for _ in range(n_nodes)]))

rhos = simulation_loop(td, states, n_nodes, 0.3)
'''