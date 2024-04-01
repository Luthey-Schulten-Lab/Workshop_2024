
import math as m
import numpy as np

import os
import csv
import pandas as pd
import importlib

import h5py


def speciesandIniCounts():
    """
    Input: None

    Return: list of species and their initial counts

    
    """
    species = []
    ini_counts = []

    species.append('P_0001')
    ini_counts.append(148)
    
    nonOricSpec = ['High_Affinity_Site','High_Affinity_Bound','High_Affinity_Site_oriC',
                   'High_Affinity_Bound_oriC','Low_Affinity_Site_1','Low_Affinity_Site_2',
                   'Low_Affinity_Bound_1','Low_Affinity_Bound_2','chromosome_C','chromosome_CC']
    nonOricSpecCounts = [16,0,1,0,0,0,0,0,1,1]
    species.extend(nonOricSpec)
    ini_counts.extend(nonOricSpecCounts)

    for i in range(30):         #loop adds 30 terms for unbound sites
        term = 'ssDNAunboundSite_'
        term = term + str(i+1)
        species.append(term)
        ini_counts.append(0)
    for k in range(30):         # add 30 more species for unbinding reactions
        unbnd = 'ssDNA_Unbinding_'
        unbnd = unbnd + str(k+1)
        species.append(unbnd)
        ini_counts.append(0)

    species.append('Initiator_C')
    ini_counts.append(0)
    species.append('Initiator_CC')
    ini_counts.append(0)
    # Fake marking species
    species.append('RepInitCheck')
    ini_counts.append(0)
    
    return species, ini_counts



def recordTraces(trajDF, lmfilename, time, restartInterval, species):
    """
    Input: 
        trajDF: data frame of trajectories of all species for one replicate
        lmfilename: the name of the LM file to be read, string
        time: time in the current simulation, integer in unit of minute
        restartInterval: restartInterval of CME simulation, integer in unit of seconds
        species: the list of all species
    return: 
        counts: The counts of all species at the end of each CME simulation, 1d array
        trajDF: updated data frame
    Description: 
        Update the data frame and return the final counts of all species of each CME simulation
    """

    f=h5py.File(lmfilename, "r")

    data=f['Simulations'][str(1).zfill(7)]['SpeciesCounts'][()].transpose()
    
    counts = data[:,-1]

    f.close()

    if time ==0:
        trajDF['Time'] = species
        timespan = np.arange(0, restartInterval+1)

        for i_loc, timemoment in enumerate(timespan):
            trajDF[timemoment] = data[:,i_loc]

    else:
        timespan = np.arange(time*restartInterval+1, (time+1)*restartInterval+1)
        for i_loc, timemoment in enumerate(timespan):
            trajDF[timemoment] = data[:,i_loc]

    
    return counts, trajDF


def volumeChange(time, doublingTime, cellVolume_init):
    """
    
    Description: Approximate time-dependent volume based on the Cell 2022's result
    The volume changes in a linear then constant way, approximately from the Figure 3E of Cell 2022 paper.
    """
    
    if time <= 60:
        cellVolume = cellVolume_init*(1+time/doublingTime)
    else:
        cellVolume = 2*cellVolume_init
    return cellVolume


def addRepInit(sim, cellVolume, counts, species, volumechange_flag):
    """
    Input:
        sim: CME simulation handle
        cellVolume: time depedent cell volume, float number
        species: the list of all species
        volumechange_flag: True or False, True to use the time dependent volume, False not to use,  Boolean type

    Return: 
        None

    Description:
        define and add species; add reactions for each CME simulation     


    """

    # Define initiation of one-time replication

    ##################################
    # Reactions constants
    NA = 6.022*(10**(23))

    if volumechange_flag:
        k_high = 7800*1000/NA/cellVolume # 0.386 molecule^-1 sec^-1
        k_low = 35*1000/NA/cellVolume  # 0.0017 molecule^-1 sec^-1
        k_on = 100*1000/NA/cellVolume # 0.005 molecule^-1 sec^-1 Number of P_0001/dnaA is 150 initially 
        k_off = 0.55 #sec^-1
    else:
        r_cell = 2.0*(10**-7) # m
        cellVolume_init = (4*np.pi/3)*1000*r_cell**3 # L
        k_high = 7800*1000/NA/cellVolume_init # 0.386 molecule^-1 sec^-1
        k_low = 35*1000/NA/cellVolume_init  # 0.0017 molecule^-1 sec^-1
        k_on = 100*1000/NA/cellVolume_init # 0.005 molecule^-1 sec^-1 Number of P_0001/dnaA is 150 initially 
        k_off = 0.55 #sec^-1


    # Reaction constants for DnaA removal when replication starts
    helicase_removal_rate = 600 #s^-1  SUPER FAST REMOVAL
    #########################################
    # Define species and counts
    sim.defineSpecies(species)

    for i in range(len(counts)):
        sim.addParticles(species = species[i], count = int(counts[i]))

    #################################
    # Add reactions
    # Step I:
    # DnaA(IV) Binds with dsDNA to unwind and open a ssDNA pocket
 
    # High_Affinity_Site is other sites nonrelevant to DNA replication and is competitive to the binding of Site_oriC
    sim.addReaction(('High_Affinity_Site', 'P_0001'),'High_Affinity_Bound',k_high)
    sim.addReaction(('High_Affinity_Site_oriC', 'P_0001'),('High_Affinity_Bound_oriC','Low_Affinity_Site_1'),k_high)
    sim.addReaction(('Low_Affinity_Site_1', 'P_0001'),('Low_Affinity_Bound_1','Low_Affinity_Site_2'),k_low)
    sim.addReaction(('Low_Affinity_Site_2', 'P_0001'),('Low_Affinity_Bound_2','ssDNAunboundSite_1'),k_low)
    # The opening of ssDNA pocket is not reversible
    
    ################################################
    # Step II: DnaA filament growth on ssDNA
    # The extension of DnaA on ssDNA is totally reversible
    # Define the species needed to represent DnaA filament growth
    # Filament Growth Reactions
    for i in range (1,30):
        sim.addReaction(('ssDNAunboundSite_' + str(i),'P_0001'), ('ssDNAunboundSite_' + str(i+1)),k_on)
        sim.addReaction(('ssDNAunboundSite_' + str(i+1)),('ssDNAunboundSite_' + str(i),'P_0001'),k_off)

    sim.addReaction(('ssDNAunboundSite_30', 'P_0001'),('ssDNA_Unbinding_30','Initiator_C','Initiator_CC','RepInitCheck'),k_on)
    
    ####################################### 
    # Step III: Separation of DnaA and ssDNA
    # Add the unbinding reactions for each of the 30 possible unbinding events in filament formation.
    # Propensities of removal reaction are rather large and usually happnes in a sudden
    for i in range (2,31):
        sim.addReaction(('ssDNA_Unbinding_' + str(i)),('ssDNA_Unbinding_' + str(i-1),'P_0001'),helicase_removal_rate)

    sim.addReaction(('ssDNA_Unbinding_1'),'P_0001',helicase_removal_rate)
    

    return None

