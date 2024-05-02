"""
User defined Python Script to use matplotlib to plot the distribution of protein and mRNA count
"""

import matplotlib.pyplot as plt
import os
import numpy as np

def plot_histogram(data, figure_path, bins, xlabel, title):
    
    fig_size = [87,87/1.618]

    fig = plt.figure(figsize=(fig_size[0]/25.4,fig_size[1]/25.4))

    plt.hist(data, bins=bins, alpha=0.7, color='limegreen')

    ax = plt.gca()

    xlabel = xlabel.replace('_','\_')
    ax.set_xlabel(r'{0}'.format(xlabel),
                  fontsize=7,
                  labelpad=1.5)

    ax.set_ylabel(r'{0}'.format('Frequency'),
                fontsize=7,
                labelpad=1.5)

    title = title.replace('_','\_')
    ax.set_title(r'{0}'.format(title),
                 fontsize=8,
                 pad=4)

    tick_length = 4.0
    tick_width = 1.5
    ax.tick_params(labelsize=5,
                    length=tick_length,
                    width=tick_width,
                    direction='in',
                    left=True,
                    right=True,
                    bottom=True,
                    top=True,
                    which='major')

    ax.tick_params(labelsize=5,
                    length=tick_length/1.5,
                    width=tick_width/1.5,
                    direction='in',
                    left=True,
                    right=False,
                    bottom=True,
                    top=False,
                    which='minor')

    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)

    mean = np.mean(data)
    plt.axvline(mean, color='red', linestyle='dashed', linewidth=1.5, label = 'Mean: {0:.3f}'.format(mean))

    median = np.median(data)
    plt.axvline(median, color='black', linestyle='dashed', linewidth=1.5, label = 'Median: {0:.3f}'.format(median))      

    min = np.min(data)
    plt.axvline(min, color='blue', linestyle='dashed', linewidth=1.5, label = 'Min: {0:.3f}'.format(min))  

    ax.legend(fontsize = 8)

    fig.savefig(figure_path, dpi = 300)

    plt.close()
    
    return None