#!/usr/bin/env python
# coding: utf-8

# # Tutorial 1.2 - Stochastic Gene Expression
# 
# Here we examine a CME model of stochastic Gene Expression.
# 
# In this model, we include the transcription and translation of gene and mRNA together with degradation of both mRNA and protein.
# 
# The model presented here can be found in the article: [Analytical distributions for stochastic gene expression](https://www.pnas.org/doi/full/10.1073/pnas.0803850105).
# 

# In[1]:


# Import Standard Python Libraries
import os
import numpy as np
import matplotlib.pyplot as plt

# Import pyLM Libraries
from pyLM import *
from pyLM.units import *
from pySTDLM import *
from pySTDLM.PostProcessing import *

# Enable plotting inline in the Jupyter notebook
get_ipython().run_line_magic('matplotlib', 'inline')


# ## Constants
# 
# Rates of first three first-order reactions come from the CME/ODE Whole Cell Model at initial conditions.
# 
# Degradation rate of protein is calculated based on 25 hours' half life in [Maier et al, 2011](https://www.embopress.org/doi/full/10.1038/msb.2011.38).

# In[2]:


# Constants
v0  = 6.41e-4       # Transcription, s^-1
d0 = 2.59e-3     # degradation of mRNA, s^-1
v1 = 7.2e-2        # translation, s^-1
d1 = 7.70e-6      # degradation of protein, s^-1


# ## Define CME simulation

# We begin by creating a [CMESimulation](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pyLM.CME.html#module-pyLM.CME) "object" that we call ```sim```. This object will include the definition of the whole stochastic simulation.

# In[3]:


# Create our CME simulation object
sim = CME.CMESimulation(name='Gene Expression')


# Next we define the chemical species with simulation. First. we specify the names of the chemical species.  Then we register these species with the simulation.  The [```defineSpecies()```](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pyLM.CME.html#module-pyLM.CME) function can be called multiple times and will add any new names to the list of species.

# In[4]:


# Define our chemical species
species = ['gene', 'mRNA', 'ptn']
sim.defineSpecies(species)


# Here we add reactions to the simulation. We use the [```addReaction()```](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pyLM.CME.html#module-pyLM.CME) function that is a member of the ```CMESimulation``` object. We add a bimolecular association reaction and a unimolecular dissociation reaction. When more than one reactant is involved, the list of reactant names should be passed as a tuple as can be seen in the reactant of the association reaction. The rates in this command must be in units of molecules and seconds, for instance units of ```/molecule/sec``` for the association reaction.

# In[5]:


# Add reactions to the simulation

sim.addReaction(reactant='gene', product=('gene','mRNA'), rate=v0)
sim.addReaction(reactant='mRNA', product='', rate=d0)
sim.addReaction(reactant='mRNA', product=('mRNA','ptn'), rate=v1)
sim.addReaction(reactant='ptn', product='', rate=d1)


# Next, we add the initial particle counts to the simulation using the [```addParticles()```](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pyLM.CME.html#module-pyLM.CME) function.

# In[6]:


# Set our initial species counts

sim.addParticles(species='gene', count=1)
sim.addParticles(species='mRNA', count=1)
sim.addParticles(species='ptn', count=0)


# Finally, we define the simulation execution parameters. We have the simulation run for 6300 seconds of real time to cover the entire cell cyle.
# 
# The traces are recorded per 1 second.
# 
# Then we name the simulation output file and save the simulation definition to it.

# In[7]:


# Simulation time is 6300, entire cell life cycle.
writeInterval = 1
simtime = 6300

sim.setWriteInterval(writeInterval)
sim.setSimulationTime(simtime)

filename = "./T2.1-geneExpression.lm"

os.system("rm -rf %s"%(filename)) # Remove previous LM file 

sim.save(filename)


# In[8]:


# Print out the information of the system
sim


# ## Run Simulation

# In[9]:


# Run multipled replicates using the Gillespie solver
# Will cost less than 1 minute
reps = 200

sim.run(filename=filename, method="lm::cme::GillespieDSolver", replicates=reps)


# ## Post-Processing

# In[10]:


fileHandle = PostProcessing.openLMFile(filename) # Create h5py file handle

plotfolder = './plots_gene_expression/'

if not os.path.exists(plotfolder):
    os.mkdir(plotfolder)


# #### Using [```plotAvgVarFromFile()```](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pySTDLM.PostProcessing.html#module-pySTDLM.PostProcessing) built-in function in PostProcessing to plot the average and variance of species. 

# In[11]:


for i_specie, specie in enumerate(species):

    picturepath = plotfolder + 'Trace.{0}.png'.format(specie)

    PostProcessing.plotAvgVarFromFile(filename = filename, species = [specie], outfile = picturepath)


# #### Using user-defined function to serialize traces in LM file to a 3D Numpy Array with dimesions *(reps, species, time)*.

# In[12]:


import h5py

timestep = PostProcessing.getTimesteps(fileHandle) # use PostProcessing to get the timesteps of the simulation

traces = np.zeros((reps,len(sim.particleMap),len(timestep)))

def get_sim_data(filename):

    f=h5py.File(filename, "r")

    for r in range(reps):
        
        traces[r]=f['Simulations'][str(r+1).zfill(7)]['SpeciesCounts'][()].transpose()
        
    f.close()

    return traces

traces = get_sim_data(filename)

# print(traces)


# #### Plot the distribution of Protein at the end of the simulation

# User defined function to draw histogram

# In[13]:


def plot_histogram(data, figure_path, bins, xlabel, title):
    
    fig_size = [87,87/1.618]

    fig = plt.figure(figsize=(fig_size[0]/25.4,fig_size[1]/25.4))

    plt.hist(ptns_end, bins=20, alpha=0.7, color='limegreen')

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

    mean = np.mean(ptns_end)
    plt.axvline(mean, color='red', linestyle='dashed', linewidth=1.5, label = 'Mean: {0:.3f}'.format(mean))

    median = np.median(ptns_end)
    plt.axvline(median, color='black', linestyle='dashed', linewidth=1.5, label = 'Median: {0:.3f}'.format(median))      

    min = np.min(ptns_end)
    plt.axvline(min, color='blue', linestyle='dashed', linewidth=1.5, label = 'Min: {0:.3f}'.format(min))  

    ax.legend(fontsize = 8)

    fig.savefig(fig_path, dpi = 300)

    plt.close()
    
    return None


# In[14]:


ptns_end = traces[:,2,-1] # Slice Numpy array to get the counts of ptn at the end of the whole cell cycle
# print(ptns_end)

fig_path = plotfolder + 'Distribution_Ptns_End.png'

xlabel = 'Counts of Protein [#]'

title = 'Distribution of Ptn Counts at the end of the simulation'
        
plot_histogram(data = ptns_end, figure_path = fig_path, bins = 20, xlabel = xlabel, title = title)

# plt.hist(ptns_end)


# In[15]:


ptns_end = traces[:,1,-1] # Slice Numpy array to get the counts of ptn at the end of the whole cell cycle
# print(ptns_end)

fig_path = plotfolder + 'Distribution_mRNA_End.png'

xlabel = 'Counts of mRNA [#]'

title = 'Distribution of mRNA Counts at the end of the simulation'
        
plot_histogram(data = ptns_end, figure_path = fig_path, bins = 5, xlabel = xlabel, title = title)


# In[16]:


break


# In[ ]:


from scipy.special import gamma

a = v0/d1
b = v1/d0

print('a',a,'b',b)
ns = np.linspace(0, 250, 251)

prefactor = gamma(ns+a)/gamma(ns+1)/gamma(a)

print(prefactor)


# In[ ]:


import math as m

m.gamma(3000)


# ## Questions 2.1
# 
# 1. Do mRNA and protein reach steady-state during the entire cell cycle? How can you tell this from the plots?
# 
# 2. In the whole cell model of JCVI-SYN3A, the initial count of protein (P_0001/DnaA) is 148. Does the protein count double during the entire cell cyle? And why this is important?
# 
# 3. (Challenge) Try to increase the simulation time to see the steady state of protein.
# 
# 

# In[ ]:




