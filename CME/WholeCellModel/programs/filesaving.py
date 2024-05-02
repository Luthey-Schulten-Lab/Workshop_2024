"""
Author: Enguang Fu

Date: March 2024

export the time traces of counts, Surface area, and fluxes into CSV files
"""

import pandas as pd
import numpy as np

import  integrate
import rxns_ODE as ODE

from math import floor
from math import log10

def writeCountstoCSV(restartNum, sim_properties):
    """

    Description:  Write the trajectories of all species into CSV file per restartCME
                Also filling the possible shorten trajectories in sim_properties['counts'] due to the jumping of large waiting time by the last value of each trajectory 
    """
    
    countsfile_path = sim_properties['path']['counts']

    restartInterval = sim_properties['restartInterval']
    
    hookInterval = sim_properties['hookInterval']

    currenttime_second = sim_properties['time_second'][-1]


    if restartNum == 0:
        particleDF = pd.DataFrame()

        spec_IDs = []

        for specID in sim_properties['counts'].keys():
            spec_IDs.append(specID)
        
        particleDF['Time'] = spec_IDs

        hookMoments = np.arange(restartNum*restartInterval,(restartNum+1)*restartInterval+hookInterval, hookInterval)

        for hookMoment in hookMoments:
            
            # location of trajectory in countsDic
            i_loc = int(hookMoment/hookInterval)

            if hookMoment > currenttime_second:
                print('Time moment {0} second is not covered between {1} second and {2} second due to a long reaction waiting time'.
                      format(hookMoment, restartNum*restartInterval,(restartNum+1)*restartInterval) )
                for count in sim_properties['counts'].values():
                    count.append(count[-1])
                print('The counts of species at time {0} second are repeated by count at time {1} second '.format(hookMoment, currenttime_second))

            new_counts = []
            for species, count in sim_properties['counts'].items():
                # 1) Certain counts are not updated during the simulation so use try-expect to avoid list index out of range
                try:
                    new_counts.append(count[i_loc])
                except:
                    new_counts.append(count[-1])

            particleDF[hookMoment] = new_counts
        particleDF.to_csv(countsfile_path,index=False)

    else:
        hookMoments = np.arange(restartNum*restartInterval+hookInterval,(restartNum+1)*restartInterval+hookInterval, hookInterval)

        particleDF = pd.read_csv(countsfile_path)

        for hookMoment in hookMoments:
    
            i_loc = int(hookMoment/hookInterval)

            if hookMoment > currenttime_second:
                print('Time moment {0} second is not covered between {1} second and {2} second due to a long reaction waiting time'.
                      format(hookMoment, restartNum*restartInterval,(restartNum+1)*restartInterval) )

                for count in sim_properties['counts'].values():
                    count.append(count[-1])
                print('The counts of species at time {0} second are repeated by count at time {1} second '.format(hookMoment, currenttime_second))

            new_counts = []
            for count in sim_properties['counts'].values():
                # Certain counts are not updated per second so use try-expect to avoid list index out of range
                try:
                    new_counts.append(count[i_loc])
                except:
                    new_counts.append(count[-1])

                
            particleDF[hookMoment] = new_counts

        particleDF.to_csv(countsfile_path,index=False)

    return None

def writeSAtoCSV(restartNum, sim_properties):
    """
    
    Description: Write the trajectories of Surface area and volume into CSV file per restartCME
        Also filling the possible shorten trajectories in sim_properties['SA'] due to the jumping of large waiting time by the last value of each trajectory 

    """

    SAfile_path = sim_properties['path']['SA']


    restartInterval = sim_properties['restartInterval']
    
    hookInterval = sim_properties['hookInterval']

    currenttime_second = sim_properties['time_second'][-1]

    if restartNum == 0:
        SADF = pd.DataFrame()

        IDs = []

        for specID in sim_properties['SA'].keys():
            IDs.append(specID)
        
        IDs.append('volume_L')

        SADF['Time'] = IDs

        hookMoments = np.arange(restartNum*restartInterval,(restartNum+1)*restartInterval+hookInterval, hookInterval)

        for hookMoment in hookMoments:
            
            i_loc = int(hookMoment/hookInterval)

            if hookMoment > currenttime_second:

                # print('Time moment {0} second is not covered between {1} second and {2} second due to a long reaction waiting time'.
                      # format(hookMoment, restartNum*restartInterval,(restartNum+1)*restartInterval) )
                #sim_properties['time_second'].append(hookMoment)
                #print('{0} second is appened to the hookMoment')
                for number in sim_properties['SA'].values():
                    number.append(number[-1])
                sim_properties['volume_L'].append(sim_properties['volume_L'][-1])

                print('The SA and volume at time {0} second are repeated by that at time {1} second'.format(hookMoment, currenttime_second))
            # 0, 1, 2, ..., 60 second
            new_numbers = []
            for numbers in sim_properties['SA'].values():
                new_numbers.append(numbers[i_loc])

            new_numbers.append(sim_properties['volume_L'][i_loc])

            SADF[hookMoment] = new_numbers

        SADF.to_csv(SAfile_path,index=False)

    else:
        particleDF = pd.read_csv(SAfile_path)

        hookMoments = np.arange(restartNum*restartInterval+hookInterval,(restartNum+1)*restartInterval+hookInterval, hookInterval)

        for hookMoment in hookMoments:
            # 61, 62, ..., 120 second, ...
            
            i_loc = int(hookMoment/hookInterval)

            if hookMoment > currenttime_second:

                for number in sim_properties['SA'].values():
                    number.append(number[-1])
                sim_properties['volume_L'].append(sim_properties['volume_L'][-1])

                print('The SA and volume at time {0} second are repeated by that at time {1} second'.format(hookMoment, currenttime_second))
            
            new_numbers = []
            for numbers in sim_properties['SA'].values():
                new_numbers.append(numbers[i_loc])

            new_numbers.append(sim_properties['volume_L'][i_loc])

            particleDF[hookMoment] = new_numbers

        particleDF.to_csv(SAfile_path,index=False)


    return None







def round_sig(x, sig=2):
    negative = False
    if x < 0:
        negative = True
    x = abs(x)
    if negative:
        return -1*round(x, sig-int(floor(log10(abs(x))))-1)
    elif x==0.0:
        return 0.0
    else:
        return round(x, sig-int(floor(log10(abs(x))))-1)


def writeFluxtoCSV(restartNum, sim_properties):
    """
    
    Called after the finish of one CME simulation

    Description: Write the fluxes through ODE reactions per restartCME
    # Since the kinetic parameters remain during the simulation, we can generate the fluxes based on one odecell object
    # The difficuly for per CMEInterval is how to get the concentration list of metabolites
    # We will create a new subdictionary called sim_properties['conc'] to store the concetrations list
    # the conc subdictionary starts from hookInterval seconds not 0 seconds
    """
    # initialize odecell again
    odemodel = ODE.initModel(sim_properties)

    Fluxfile_path = sim_properties['path']['flux']

    concDict = sim_properties['conc']
    # print(concDict)
    solver = integrate.noCythonSetSolver(odemodel)

    restartInterval = sim_properties['restartInterval']
    
    hookInterval = sim_properties['hookInterval']

    currenttime_second = sim_properties['time_second'][-1]

    hookMoments = np.arange(restartNum*restartInterval+hookInterval,(restartNum+1)*restartInterval+hookInterval, hookInterval)

    if restartNum == 0:
        FluxDF = pd.DataFrame()

        rxnIDs = []
        for index,rxn in enumerate(odemodel.getRxnList()):
            rxnIDs.append('F_' + rxn.getID())

        FluxDF['Time'] = rxnIDs

        for hookMoment in hookMoments:
            # 1, 2, 3, ... restartInterval
            
            i_loc = int(hookMoment/hookInterval)

            if hookMoment > currenttime_second:
                for conc in concDict.values():

                    conc.append(conc[-1])
                print('The concentration of metabolites at time {0} second are repeated by that at time {1} second'.format(hookMoment, currenttime_second))
            
            currentConc = []
            for conc in concDict.values():
                currentConc.append(conc[i_loc-1])

            currentFlux = solver.calcFlux(0, currentConc)
            new_flux= []

            for index, rxn in enumerate(odemodel.getRxnList()):

                new_flux.append((round_sig(currentFlux[index], sig=3)))

            FluxDF[hookMoment] = new_flux

        FluxDF.to_csv(Fluxfile_path,index=False)

    else:
        FluxDF = pd.read_csv(Fluxfile_path)

        for hookMoment in hookMoments:

            i_loc = int(hookMoment/hookInterval)

            if hookMoment > currenttime_second:
                for conc in concDict.values():

                    conc.append(conc[-1])
                print('The concentration of metabolites at time {0} second are repeated by that at time {1} second'.format(hookMoment, currenttime_second))

            currentConc = []
            for conc in concDict.values():
                currentConc.append(conc[i_loc-1])

            currentFlux = solver.calcFlux(0, currentConc)
            new_flux = []
            
            for index, rxn in enumerate(odemodel.getRxnList()):

                new_flux.append((round_sig(currentFlux[index], sig=3)))

            FluxDF[hookMoment] = new_flux

        FluxDF.to_csv(Fluxfile_path,index=False)

    return None

def appendHookMoments(restartNum, sim_properties):
    """
    Description:
        To make time moments in sim_properties complete by appending time points
        solve the issue that the possible jumping of CME time over several hookIntervals 
    """

    restartInterval = sim_properties['restartInterval']
    
    hookInterval = sim_properties['hookInterval']

    hookMoments = np.arange(restartNum*restartInterval,(restartNum+1)*restartInterval+hookInterval, hookInterval)

    currenttime_second = sim_properties['time_second'][-1]

    for hookMoment in hookMoments:

        if hookMoment > currenttime_second:

            sim_properties['time_second'].append(hookMoment)
            
            print('{0} second is appened to the simulationTime'.format(hookMoment))

    return None

