"""
Author: Enguang Fu

Date: March 2024

the operations when do hooking/communicating between CME/ODE
"""

import time as phys_time


import communicate
import integrate
import rxns_ODE as ODE
import filesaving
import initiation as IC


def hook_CMEODE(sim_properties, genome3A):
    """
    Input: sim_properties, genome3A
    
    Description: 1) calculate the cost of nucleotides and pass the extracted counts of ODE species to odecell;
                 2) Do ODE simulation
                 3) Record the ODE counts into counts and conc dictionaries in sim_properties and update ODE result to CME
                 3) update the Surface area
    """
    
    
    # calculate the cost of gene expression
    communicate.calculateCosts(sim_properties, genome3A, sim_properties['LocusNumtoIndex'])
    # update the counts of CME species that are also in ODE
    communicate.communicateCostsToMetabolism(sim_properties)


    iniode_start = phys_time.time()

    # Initialize/construct the ODEs using odecell
    # The counts in the counts dictionary are passed to ODE solver
    odemodel = ODE.initModel(sim_properties)

    print('ODE object Initialized in {0} seconds'.format(phys_time.time() - iniode_start))

    runode_start = phys_time.time()

    initVals = integrate.getInitVals(odemodel)

    solver = integrate.noCythonSetSolver(odemodel)

    odelength = sim_properties['hookInterval']

    odeResults = integrate.runODE(initVals, solver, odemodel, odelength)

    print('ODE simulation of {0} second Finished in {1} seconds'.format(odelength, phys_time.time()- runode_start))

    # Create conc dictionary only for the calculation of flux
    if sim_properties['time_second'][-1] == 0:
        IC.initializeConcDictionary(sim_properties, odemodel)
    
    # Save the concentrations directly from ODE to conc dictionary 
    communicate.saveConc(sim_properties, odeResults, odemodel)
    
    # update the counts of ODE species into counts dictionary
    communicate.updateODEcounts(sim_properties, odeResults, odemodel)

    # update the Surface Area and volume
    communicate.updateSA(sim_properties)    
    
    sim_properties['time_second'].append(sim_properties['time_second'][-1] + sim_properties['hookInterval'])


    return None