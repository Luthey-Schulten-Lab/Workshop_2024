"""
Author: Enguang Fu

Date: March 2024

user defined solver where hookSimulation performs every hookInterval
"""

import communicate

import hook_CMEODE

# import pyLM.units to make sure lm.GillespieDSolver is readabe by python
from pyLM.units import *

class MyOwnSolver(lm.GillespieDSolver):

    def initializeSolver(self, sim, CMECounts, sim_properties,genome3A):
        self.sim = sim
        self.CMECounts = CMECounts
        self.sim_properties = sim_properties
        self.genome3A = genome3A



    def hookSimulation(self, time):
        
        restartNum = self.sim_properties['restartNum'][-1]
        restartInterval = self.sim_properties['restartInterval']
      
        print('*******************************************************************')
        print('HookSimulation is called at '+ str(restartNum*restartInterval + time) + ' second')
        
        self.CMECounts.update(self)
        # using self.CMECounts to update CME counts and write them into sim_properties which will be used in ODE
        communicate.updateCMEcountsHook(self.sim, self.CMECounts, self.sim_properties)
        
        # Communicate Cost; Do ODE; Update SA and volume
        hook_CMEODE.hook_CMEODE(self.sim_properties, self.genome3A)
        
        # update counts of certain ODE species that also in CME to self.CMECounts
        communicate.updateODEtoCME(self.sim, self.CMECounts, self.sim_properties)

        return 1
        


        


            