# Import Standard Python Libraries
import os
import numpy as np
import matplotlib.pyplot as plt

# Import pyLM Libraries
from pyLM import *
from pyLM.units import *
from pySTDLM import *
from pySTDLM.PostProcessing import *


# Constants
v0  = 6.41e-4       # Transcription, s^-1
d0 = 2.59e-3     # degradation of mRNA, s^-1
v1 = 7.2e-2        # translation, s^-1
d1 = 7.70e-6      # degradation of protein, s^-1

# Create our CME simulation object
sim = CME.CMESimulation(name='Gene Expression')

# Define our chemical species
species = ['gene', 'mRNA', 'ptn']
sim.defineSpecies(species)

# Add reactions to the simulation
sim.addReaction(reactant='gene', product=('gene','mRNA'), rate=v0)
sim.addReaction(reactant='mRNA', product='', rate=d0)
sim.addReaction(reactant='mRNA', product=('mRNA','ptn'), rate=v1)
sim.addReaction(reactant='ptn', product='', rate=d1)

# Set our initial species counts
sim.addParticles(species='gene', count=1)
sim.addParticles(species='mRNA', count=1)
sim.addParticles(species='ptn', count=0)

# Simulation time is 6300, entire cell life cycle.
writeInterval = 1
simtime = 600

sim.setWriteInterval(writeInterval)
sim.setSimulationTime(simtime)

filename = "./T2.1-geneExpression.lm"

os.system("rm -rf %s"%(filename)) # Remove previous LM file 

sim.save(filename)

# Run multipled replicates using the Gillespie solver
# Will cost less than 1 minute
reps = 10


sim.run(filename=filename, method="lm::cme::GillespieDSolver", replicates=reps)

# Post-processing
fileHandle = PostProcessing.openLMFile(filename) # Create h5py file handle

plotfolder = './plots_gene_expression/'

if not os.path.exists(plotfolder):
    os.mkdir(plotfolder)

# Plot the average and variance of mRNA and protein
for i_specie, specie in enumerate(species):

    if specie != 'gene':
        picturepath = plotfolder + 'Trace.{0}.png'.format(specie)
        PostProcessing.plotAvgVarFromFile(filename = filename, species = [specie], outfile = picturepath)

