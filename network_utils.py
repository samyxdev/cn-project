import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import poisson
import math

# Plotting auxiliary functions

# Plot degrees probability distributions in log-log scale
def plot_distrib_log(graph, colour='#40a6d1', alpha=.8, fit_line=False, gam = 3,  expct_lo=1, expct_hi=10, expct_const=1):
    num_nodes = graph.number_of_nodes()

    # Calculate the maximum degree to know the range of x-axis
    max_degree = np.max([val for key,val in graph.degree()])

    # X-axis and y-axis values
    x = []
    y_tmp = []

    # Loop over all degrees until the maximum to compute the portion of nodes for that degree
    for i in range(max_degree + 1):
        x.append(i) # build x list
        y_tmp.append(0) # build y list

        for n in graph.nodes():
            if graph.degree(n) == i:
                y_tmp[i] += 1
        y = [i / num_nodes for i in y_tmp]

    plt.xscale('log')
    plt.yscale('log')
    plt.title('Degree distribution (log-log scale)')
    plt.ylabel('log(P(k))')
    plt.xlabel('log(k)')
    plt.plot(x, y, linewidth = 0, marker = 'o', markersize = 8, color = colour, alpha = alpha)

    if fit_line:
        # Add theoretical distribution line k^-3
        # Note that you need to parametrize it manually
        w = [a for a in range(expct_lo,expct_hi)]
        z = []
        for i in w:
            x = (i**-gam) * expct_const # set line's length and fit intercept
            z.append(x)

        plt.plot(w, z, 'k-', color='#7f7f7f')

    plt.legend(["Real degrees", f"Powerlaw model (gamma= {-gam})"])


# Plot degrees probability distributions in linear scale
def plot_distrib_lin(graph, colour='#40a6d1', alpha=.8, fit_poisson=False, fit_lambda=None, fit_binom_p=None, expct_lo=1, expct_hi=10):
    num_nodes = graph.number_of_nodes()

    # Calculate the maximum degree to know the range of x-axis
    max_degree = np.max([val for key,val in graph.degree()])

    # X-axis and y-axis values
    x = []
    y_tmp = []

    # Loop over all degrees until the maximum to compute the portion of nodes for that degree
    for i in range(max_degree + 1):
        x.append(i) # build x list
        y_tmp.append(0) # build y list

        for n in graph.nodes():
            if graph.degree(n) == i:
                y_tmp[i] += 1
        y = [i / num_nodes for i in y_tmp]

    plt.plot(x, y, linewidth = 0, marker = 'o', markersize = 8, color = colour, alpha = alpha)
    plt.title('Degree distribution (linear scale)')
    plt.ylabel('P(k)')
    plt.xlabel('k')

    if fit_poisson:
        # Add theoretical distribution poisson line
        w = [a for a in np.arange(expct_lo,expct_hi)]
        z = []

        #for i in w:
        #   x = ((lamb**i) * math.exp(-lamb)) / math.factorial(i) # poisson pmf
        #    z.append(x)

        #create an array with Poisson probability values
        z = poisson.pmf(w, mu=fit_lambda)
        plt.plot(w, z, 'k-', color='#7f7f7f')

        plt.legend(["Real degrees", f"Poisson model (lambda={fit_lambda})"])

    elif fit_binom_p != None:
        # Add theoretical distribution binomial line
        w = [a for a in np.arange(expct_lo,expct_hi)]
        z = [math.comb(num_nodes-1, a) * fit_binom_p**a * (1-fit_binom_p)*(num_nodes-1-a) for a in w]

        print(z)

        #create an array with Binomial probability values
        plt.plot(w, z, 'k-', color='#7f7f7f')

        plt.legend(["Real degrees", f"Binomial model (lambda={fit_lambda})"])

    else:
        plt.legend(["Real degrees"])


"""
        LINEAR AND LOGLOG PLOTS OF PDF AND CCDF
"""

# Auxiliary function to create logaritmically spaced bins (for log-log histogram) by specifying the number of bins
def create_log_bins(degrees, num = 20):
    bins = np.logspace(np.log10(np.min(degrees)), np.log10(np.max(degrees)), num)
    bins = np.array(bins)
    return bins

# PDF histogram in linear scale
def plot_linear_PDF(G, name='', nbins = 15):
    degrees = [G.degree(n) for n in G.nodes()]
    #plt.figure(figsize=(12,8))
    plt.title('PDF in linear scale', fontsize=15)
    plt.xlabel('Degree', fontsize=13)
    plt.ylabel('PDF', fontsize=13)
    plt.hist(degrees, bins=nbins, density = True, cumulative = False)
    plt.tight_layout()
    plt.style.use('seaborn-darkgrid')

# PDF histogram in Log-Log scale
def plot_loglog_PDF(G, name="", nbins=20):
    degrees = [G.degree(n) for n in G.nodes()]

    # creating logaritmically spaced bins
    bins = create_log_bins(degrees, num = nbins)

    #plt.figure(figsize=(12,8))
    plt.title('PDF in log-log scale',  fontsize=15)
    plt.xlabel('Degree', fontsize=13)
    plt.ylabel('PDF', fontsize=13)
    plt.yscale('log')
    plt.xscale('log')
    plt.hist(degrees, bins=bins, density = True, cumulative = False)
    plt.tight_layout()
    plt.style.use('seaborn-darkgrid')

# CCDF histogram in linear scale
def plot_linear_CCDF(G, name="", nbins=30):
    degrees = [G.degree(n) for n in G.nodes()]
    #plt.figure(figsize=(12,8))
    plt.title('CCDF in linear scale', fontsize=15)
    plt.xlabel('Degree', fontsize=13)
    plt.ylabel('CCDF', fontsize=13)
    plt.hist(degrees, bins=nbins, density = True, cumulative = -1)
    plt.tight_layout()
    plt.style.use('seaborn-darkgrid')

# CCDF histogram in Log-Log scale
def plot_loglog_CCDF(G, name="", nbins=30):
    degrees = [G.degree(n) for n in G.nodes()]

    # creating logaritmically spaced bins
    bins = create_log_bins(degrees, num=nbins)

    #plt.figure(figsize=(12,8))
    plt.title('CCDF in log-log scale', fontsize=15)
    plt.xlabel('Degree', fontsize=13)
    plt.ylabel('CCDF', fontsize=13)
    plt.yscale('log')
    plt.xscale('log')
    plt.hist(degrees, bins=bins, density = True, cumulative = -1)
    plt.tight_layout()
    plt.style.use('seaborn-darkgrid')