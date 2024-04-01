# pyLM to construct the CME model
from pyLM import *
from pyLM.units import *


import math as m
import numpy as np
import sys
import os
import csv
import pandas as pd
import importlib

import h5py

import func_initiation as f

log = open('./initiation_CME.log', 'w')

original_stdout = sys.stdout
original_stderr = sys.stderr

sys.stdout = log
sys.stderr = log

#In each simulation of whole cell cycle, restart the simulation per 60 seconds to update the cell volume
restartInterval = 60 # s

# Simulation 50 replicates serially
reps = 50

# Initial value of cell radius is 200 nm
r_cell = 2.0*(10**-7) # m
cellVolume_init = (4*np.pi/3)*1000*r_cell**3 # L

# Time needed to double the cell volume is 60 minutes
doublingTime = 60 # min

# The length of whole cell cycle is 105 minutes.
endtime = 105  # minutes

# Two-layer loop to simulate 20 replicates's 105 minute cell cycle
csvfilefolers = ['./WithVolumeChange/', './WithoutVolumeChange/']
volumeChange_flag = [True, False]

for i_change, change_flag in enumerate(volumeChange_flag):

    csvfilefolder = csvfilefolers[i_change]

    if not os.path.exists(csvfilefolder):
        os.makedirs(csvfilefolder)

    # Get the species list and initial counts
    species, ini_counts = f.speciesandIniCounts()

    for rep in range(1,reps+1):
        print('Start the simulation of '+str(rep) +' replicate')
        
        counts = ini_counts
        
        trajDF = pd.DataFrame()
        
        for time in range(0,endtime):

            # Filename of the current simulation
            filename = 'Ini_' + str(time+1) + '.lm'
            # Filename of the last minute simulation
            filename_lastmiute = 'Ini_' + str(time) +'.lm'

            # Linear increase to double volume then contains
            cellVolume = f.volumeChange(time, doublingTime, cellVolume_init)
            
            # Initialize each minute's simulation
            sim=CME.CMESimulation(name="DNA Replication Initiation in Syn3A" + str(time+1))
            # Write the trajectories of each species per 1 second
            sim.setWriteInterval(1)
            sim.setSimulationTime(restartInterval)
            
            # define the reactions and particles numbers every simulation
            f.addRepInit(sim, cellVolume, counts, species, change_flag)

            # Remove the lm file of last minute so that you will not get 105 files 
            os.system("rm -rf %s"%(filename_lastmiute))

            # Remove lm file from last whole cell cycle simulation
            os.system("rm -rf %s"%(filename))
        
            sim.save(filename)
            sim.run(filename, "lm::cme::GillespieDSolver", 1)
            
            # counts is the number of each species at last frame written time
            # trajDF is the record of numbers of all frame written time in all simulation
            counts,trajDF = f.recordTraces(trajDF, filename, time, restartInterval, species)

            #print('Simulation Time Finished ' + str(time+1) + ' Min')
        print('Finish the simulation of '+str(rep) +' replicate')
        
        # Write the trajectory to csv file
        csvpath = '{0}initiation_{1}.csv'.format(csvfilefolder, rep)
        trajDF.to_csv(csvpath, index = False)


sys.stdout = original_stdout 
sys.stderr = original_stderr 