import numpy as np
import mplhep as hep
import matplotlib.pyplot as plt
plt.style.use(hep.style.ROOT)


def plot_overall_interference(terms, names,
                              output_directory, output_filename):
    """This function plots the overall plot of both interference and pure terms to plot everything

    Arguments:
        terms -- A list of all the pure sample and interference terms. This should be a list of (count, bin) pairs (i.e. numpy histograms)
        names -- A list of the names for all of these terms
        output_directory -- The directory you would like to output to
        output_filename -- The filename you want to name the plots

    Returns:
        the numpy histogram object of the overall sample
    """
    
    if output_directory[-1] != '/':
        output_directory += '/'
    
    plt.cla()
    plt.gca().axhline(zorder=-1, color='black', lw=2)
    
    bins = terms[0][1]
    overall = np.zeros(len(bins) - 1, dtype=float)
    for n, term in enumerate(terms):
        overall += term[0]
        hep.histplot(term, label=names[n], lw=2)
        
    hep.histplot(overall, bins, label="Overall", lw=3)
    plt.legend(loc="upper right")
    plt.savefig(output_directory + output_filename + "overall_interference_plot.png")
    plt.savefig(output_directory + output_filename + "overall_interference_plot.pdf")
    
    return overall, bins