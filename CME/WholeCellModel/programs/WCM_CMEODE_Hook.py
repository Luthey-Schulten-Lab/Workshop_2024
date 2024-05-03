"""
Author: Enguang Fu
Date: March 2024

Main Scirpt for CMEODE Hybrid Whole Cell Simulation
"""

###################################
#####   Import ####
###################################

# pyLM to construct the CME system
from pyLM import *

# 
import numpy as np
import os, sys
import time as phys_time
from Bio import SeqIO

from datetime import datetime
import argparse

from mpi4py import MPI

# User defined Python scripts
import initiation as IC
import rxns_CME
import communicate
import species_counts
import hookSolver_CMEODE
import filesaving
import hook_CMEODE



###################################
#####   Parsing Parameters ####
###################################


comm = MPI.COMM_WORLD
rank = comm.Get_rank()+1
size = comm.Get_size()


# Arguments needed: simulationType, simulation time length, restartInterval, hookInterval, writeInterval, outputfolder
ap = argparse.ArgumentParser()


ap.add_argument('-st','--simType',required=True)
ap.add_argument('-t','--simTime',required=True)
ap.add_argument('-rs','--restartInterval',required=True)
ap.add_argument('-hi', '--hookInterval', required = True)
ap.add_argument('-f', '--folder', required= True)

args = ap.parse_args()

dir_path = args.folder

# Check if the directory exists
if not os.path.exists(dir_path):
        # The directory does not exist, create it
        os.makedirs(dir_path)


logfile = 'log_{0}.txt'.format(rank)

logfile_path = os.path.join(dir_path, logfile)

# Redirect the output and possible error to a log file
log = open(logfile_path,'w')

original_stdout = sys.stdout
original_stderr = sys.stderr

sys.stdout = log
sys.stderr = log


print('The PID is ' + str(os.getpid()))

try:
    
    simTime = float(args.simTime)
    restartInterval = float(args.restartInterval)
    hookInterval = float(args.hookInterval)

except:
    
    sys.exit("Error: Please enter simulation time (-t), restart time (-rs), hook interval (-hi), and write interval (-wi) as the valued time in seconds.")

if simTime % restartInterval == 0:
        if restartInterval % hookInterval == 0:
                    if (str(args.simType) == "cme-ode"):
                           None
                    else:
                                sys.exit("Error: Enter 'cme-ode' for a hybrid CME-ODE simulation")
        else:
                sys.exit('Error: Please enter restart time (-rs) as integer multiples of hook interval (-hi)')
else:
        sys.exit("Error: Please enter simulation time (-t) as integer multiples of restart time (-rs)")


restartNums = int(simTime/restartInterval)

print('The simulation time is {0} seconds with restart interval {1} seconds.'.format(simTime, restartInterval))
print('The hook interval is {0} seconds and write interval to the LM file {1} seconds'.format(hookInterval, hookInterval))
print('CME simulation will restart {0} times to finish the simulation time'.format(restartNums))


print('*******************************************************************************')


###################################
#####   Initialize sim_properties ####
###################################

start_time = datetime.now()

print('The simulation starts at '+ str(start_time))

# Initialize sim_properties, a dictionary record the genome information, trajectories, and other constants

sim_properties = {}

# Input Gene bank file syn3A.gb
genomeFile3A =  '../input_data/syn3A.gb'
genome3A = next(SeqIO.parse(genomeFile3A, "gb"))

# Convert Gene bank file into dictionary genome
IC.initializeConstants(sim_properties, genome3A)

# Initialize the counts and concentrations of metabolites
IC.initializeMetabolitesCounts(sim_properties)
IC.initializeMediumConcs(sim_properties)
IC.initializeProteinMetabolitesCounts(sim_properties)

# Initialize the counts of cost marker
IC.initializeCosts(sim_properties)


# The names of three csv files e.g. counts_1.csv, SA_2.csv, Flux_20.csv
countsCSV = 'counts_{0}.csv'.format(rank)
countsCSV_path = os.path.join(dir_path,countsCSV)

SACSV = 'SA_{0}.csv'.format(rank)
SACSV_path = os.path.join(dir_path, SACSV)

fluxCSV = 'Flux_{0}.csv'.format(rank)
fluxCSV_path = os.path.join(dir_path, fluxCSV)

sim_properties['path'] = {'counts': countsCSV_path, 'SA': SACSV_path, 'flux':fluxCSV_path}



# sim_properties['time_second] will append the hook moments until hits the simulation time length
sim_properties['time_second'] = [int(0)]

# simulaion parameters
sim_properties['restartInterval'] = restartInterval

sim_properties['hookInterval'] = hookInterval

sim_properties['restartNum'] = []



###################################
#####   Perform the Simulation ####
###################################


# Iterating Loop to restart CME restartNums times to finish the simulation time
for restartNum in range(0,restartNums):

    sim_properties['restartNum'].append(restartNum)
    print('*******************************************************************')

    print('Start the simulation between {0} to {1} seconds'.format(restartNum*restartInterval, (restartNum+1)*restartInterval))
    
    sim=CME.CMESimulation(name="CMEODE_" + str(restartNum+1))

    # Set the parameters of the CME
    sim = IC.initializeCME(sim, restartNum, sim_properties)
    
    
    # define the reactions and particles numbers every newly start simulation
    ini_start = phys_time.time()

    IC.addEnzymesCounts(sim, sim_properties)

    rxns_CME.addGeneticInformationReactions(sim, sim_properties, genome3A)

    IC.addGeneticInformationSpeciesCounts(sim,sim_properties)

    print('Initiation of CME reactions and counts in {0} seconds'.format(phys_time.time() - ini_start))

    # CMECounts need to be reinitialized per new CME simulation
    CMECounts = species_counts.SpeciesCounts(sim)
    hook = hookSolver_CMEODE.MyOwnSolver()

    hook.initializeSolver(sim, CMECounts, sim_properties, genome3A)

    if restartNum == 0:
        IC.initializeMembrane(sim_properties)
    
    # LMfilename of the current simulation
    LMfilename = os.path.join(dir_path, 'CME_ODE_{0}_{1}.lm'.format(rank, restartNum))

    # Remove lm file from previous simulation
    os.system("rm -rf %s"%(LMfilename))

    save_start = phys_time.time()

    sim.save(LMfilename)
    print('LM file Saved in {0} seconds'.format(phys_time.time() - save_start))

    runCMEODE_start = phys_time.time()
    
    # run CME simulation with hook of one single replicate
    sim.runSolver(filename = LMfilename,solver = hook, replicates = 1)

    # Now the CME simulation ends, we use updateCMEcountsFile to read in the counts and pass them to ODE to do the last second's ODE simulation

    print('HookSimulation is called at {0} second'.format((restartInterval*(restartNum+1))))

    communicate.updateCMEcountsFile(sim, sim_properties, LMfilename)
    
    hook_CMEODE.hook_CMEODE(sim_properties, genome3A)
    
    CSVfilesaving_start = phys_time.time()


    filesaving.writeCountstoCSV(restartNum,sim_properties)
    filesaving.writeSAtoCSV(restartNum, sim_properties)
    filesaving.writeFluxtoCSV(restartNum, sim_properties)
    filesaving.appendHookMoments(restartNum, sim_properties)
    

    print('CSV Files update per restartInterval of Counts, SA and Flux in {0} seconds'.format(phys_time.time() - CSVfilesaving_start))

    print('CMEODE run finished in {0} seconds'.format(phys_time.time() - runCMEODE_start))

    # Remove the just finished lm file
    os.system("rm -rf %s"%(LMfilename))

    print('Finish the simulation between {0} to {1} seconds'.format(restartNum*restartInterval, (restartNum+1)*restartInterval))
    print('*******************************************************************')



###################################
#####   Close the Simulation ####
###################################


end_time = datetime.now()
print('The simulation ends at {0}'.format(end_time))

print('Time (hour:minutes:seconds) taken to finish the simulation: {0}'.format(end_time - start_time))

sys.stdout = original_stdout
sys.stderr = original_stderr

print('The simulation ends at {0}'.format(end_time))

print('Time taken to finish the simulation: {0}'.format(end_time - start_time))

