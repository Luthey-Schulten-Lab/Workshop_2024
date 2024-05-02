"""
Author: Enguang Fu

Date: March 2024

updates CME and ODE states, calculate costs and update membrane
"""


import numpy as np
import h5py
from initiation import mMtoPart



def updateCMEcountsFile(sim, sim_properties, filename):
    """
    Input: sim_properties dictionary; filename of lm file

    Return None
      
    Called immediately after CME run; Works in restart CME per second

    Description: Read in the finished LM file using h5py; 
    Append the newest counts of CME species to CountsDic
    """


    f=h5py.File(filename, "r")
    # Take the last second's counts
    data=f['Simulations'][str(1).zfill(7)]['SpeciesCounts'][()].transpose()[:,-1]
    
    # PP.closeLMFile(f)
    
    f.close()
    
    countsDic = sim_properties['counts']

    for i in range(len(data)):
        speciesName = list(sim.particleMap.keys())[i]
        countsDic[speciesName].append(int(data[i]))


    return None

def updateCMEcountsHook(sim, CMECounts, sim_properties):

    """
    
    Description: Using David's IET hook pass the latest counts of CME species to dictionary
    """

    countsDic = sim_properties['counts']

    for CMEName in sim.particleMap.keys():
        countsDic[CMEName].append(CMECounts[CMEName])

    return None

def updateODEcounts(sim_properties, odeResults, odemodel):
    """
    Input: 
        sim_properties: Dictionary contains the counts trajectories
        odeResults: The concentration trajectories of metabolites
        odemodel: Constructed ode model by odecell  

    Return:
        None

    Called by:
        hookSimulation
        
    Description: Update the counts of metabolites after ODE simulation into the sim_properties dictionary
    """
    

    resFinal = odeResults[-1,:]
    
    metIDDict = odemodel.getMetDict()
    # {'M_ACP_c': 0, 'M_ACP_R_c': 1, 'M_apoACP_c': 2, ...}

    # print(metIDDict)

    countsDic = sim_properties['counts']

    # Species in appendedList have one more element compared to other metaboites because they have been appened in updataCME or communicateCoststoMetabolism
    # so different way to update the counts after the ODE run.
    
    appendedList = ['M_ppi_c', 'M_pi_c', 'M_atp_c', 'M_ctp_c', 'M_utp_c', 'M_gtp_c', 'M_adp_c', 'M_cdp_c', 'M_udp_c'
        , 'M_gdp_c', 'M_datp_c', 'M_dttp_c', 'M_dctp_c', 'M_dgtp_c', 'M_amp_c', 'M_ump_c', 'M_cmp_c', 'M_gmp_c']
    aa_list = sim_properties['aa_list']
    appendedList.extend(aa_list)

    for metID, Num in metIDDict.items():
        if not metID.endswith('_e'):
            if metID in appendedList:
                countsDic[metID][-1] = mMtoPart(resFinal[Num], sim_properties)
            
            else:
                countsDic[metID].append(mMtoPart(resFinal[Num], sim_properties))


    return None 

def updateODEtoCME(sim, CMECounts, sim_properties):
    """
    Input:

    Return: None

    Called by hookSimulation
    
    Description: After the ODE run, update certain CME species's counts (tRNA charging related species including 20 aas, atp, amp, and ppi) to the following CME simulation
    """

    # tRNA_ODE_list is the list of species that both in CME and ODE
    tRNA_ODE_list =  ['M_atp_c', 'M_amp_c', 'M_ppi_c']
    aa_list = sim_properties['aa_list']
    tRNA_ODE_list.extend(aa_list)
    countsDic = sim_properties['counts']

    for species in tRNA_ODE_list:
        CMECounts[species] = countsDic[species][-1]


    return None


def calculateCosts(sim_properties, gbfile, LocusNumtoIndex):
    """
    Input: None

    Return: None

    Called per communication step

    Description: calculate the costs of replication, transcription, translaion, degradation and translocation
    
    """

    calculateReplicationCosts(sim_properties, gbfile, LocusNumtoIndex)

    calculateTranscriptionCosts(sim_properties)

    calculateDegradationCosts(sim_properties)

    calculateTranslationCosts(sim_properties)

    calculateTranslocationCosts(sim_properties)

    calculatetRNAchargingCosts(sim_properties)

    return None


def calculateReplicationCosts(sim_properties, gbfile, LocusNumtoIndex):
    """
    Input: sim_properties dictionary; genebank file; LocusNumtoIndex dictionary

    Output: None

    Called by hookSimulation

    Description: Calcualte the dNTPs and ATP costs in the 1 second CME simulation and append the costs to corresponding cost species in sim_properties  
    
    1 ATP hydrolysis per bp to unwind the dsDNA
    """

    # Calculate the dNTP consumption based on the RNAsequence
    # The intergenic region not included; should be improved
    dna = gbfile

    genome = sim_properties['genome']

    countsDic = sim_properties['counts']

    nuc_repcost = {'dATP_DNArep_cost':'M_datp_c', 'dTTP_DNArep_cost':'M_dttp_c', 'dCTP_DNArep_cost':'M_dctp_c', 'dGTP_DNArep_cost':'M_dgtp_c'}
    # pseudoGenes = ['JCVISYN3A_0051', 'JCVISYN3A_0546','JCVISYN3A_0602']
    pseudoGenes = sim_properties['pseudoGenes']
    ATP_replication_cost = 0
    dATP_rep_cost = 0
    dTTP_rep_cost = 0
    dCTP_rep_cost = 0
    dGTP_rep_cost = 0

    for locusTag, locusDic in genome.items():
        # Three pseudo Genes are not replicated
        if locusTag not in pseudoGenes:
            locusNum = locusTag.split('_')[1]

            Produced_gene = 'Produced_G_' +locusNum
            geneGenerated = int(countsDic[Produced_gene][-1] - countsDic[Produced_gene][-2])

                # read the index of the start and end point of the genic and intergenic region from prebuilt dictionary
            index = LocusNumtoIndex[locusNum]
            # print(index)
            dnasequence = str(dna.seq[index[0]:index[1]])
                
            dATP_count = dnasequence.count('A')
            dTTP_count = dnasequence.count('T')
            dCTP_count = dnasequence.count('C')
            dGTP_count = dnasequence.count('G')
            
            # The dNTPs as monomers are used on both leading and lagging strand
            dATP_rep_cost += (dATP_count+dTTP_count)*geneGenerated
            dTTP_rep_cost += (dTTP_count+dATP_count)*geneGenerated
            dCTP_rep_cost += (dCTP_count+dGTP_count)*geneGenerated
            dGTP_rep_cost += (dGTP_count+dCTP_count)*geneGenerated

            ATP_replication_cost += (dATP_count + dTTP_count + dCTP_count + dGTP_count)*geneGenerated



    countsDic['dATP_DNArep_cost'].append(dATP_rep_cost)
    countsDic['dTTP_DNArep_cost'].append(dTTP_rep_cost)
    countsDic['dCTP_DNArep_cost'].append(dCTP_rep_cost)
    countsDic['dGTP_DNArep_cost'].append(dGTP_rep_cost)
    countsDic['ATP_DNArep_cost'].append(ATP_replication_cost)

    return None

def calculateTranscriptionCosts(sim_properties):
    """
    Input: sim_properties dictionary
    
    Return: None

    Called by hookSimulation

    Description: calculate the ATP and NTPs costs in the transcription for mRNA, tRNA and rRNA
    and append the costs to corresponding species

    1 ATP hydrolysis per bp to unwind the dsDNA
    
    """
    countsDic = sim_properties['counts']
    genome = sim_properties['genome']
    ATP_trsc_energy_cost = 0 
    ATP_mRNA_cost = 0; ATP_rRNA_cost = 0; ATP_tRNA_cost = 0
    UTP_mRNA_cost = 0; UTP_rRNA_cost = 0; UTP_tRNA_cost = 0
    CTP_mRNA_cost = 0; CTP_rRNA_cost = 0; CTP_tRNA_cost = 0
    GTP_mRNA_cost = 0; GTP_rRNA_cost = 0; GTP_tRNA_cost = 0

    nuc_trsccosts = {'ATP_mRNA_cost':'M_atp_c', 'CTP_mRNA_cost':'M_ctp_c', 'UTP_mRNA_cost':'M_utp_c', 'GTP_mRNA_cost':'M_gtp_c'}

    for locusTag, locusDic in genome.items():
        if locusDic['Type'] == 'protein':
            locusNum = locusTag.split('_')[1]
            Produced_RNA = 'Produced_R_' +locusNum
            RNAGenerated = int(countsDic[Produced_RNA][-1] - countsDic[Produced_RNA][-2])


            rnasequence = locusDic['RNAsequence']
            ATP_count = rnasequence.count('A')
            UTP_count = rnasequence.count('U')
            CTP_count = rnasequence.count('C')
            GTP_count = rnasequence.count('G')

            ATP_mRNA_cost += ATP_count*RNAGenerated
            UTP_mRNA_cost += UTP_count*RNAGenerated
            CTP_mRNA_cost += CTP_count*RNAGenerated
            GTP_mRNA_cost += GTP_count*RNAGenerated

            ATP_trsc_energy_cost += (ATP_count + UTP_count + GTP_count + CTP_count)*RNAGenerated

        elif locusDic['Type'] == 'rRNA':
            Produced_RNA = 'Produced_R_' +locusNum
            RNAGenerated = int(countsDic[Produced_RNA][-1] - countsDic[Produced_RNA][-2])


            rnasequence = locusDic['RNAsequence']
            ATP_count = rnasequence.count('A')
            UTP_count = rnasequence.count('U')
            CTP_count = rnasequence.count('C')
            GTP_count = rnasequence.count('G')

            ATP_rRNA_cost += ATP_count*RNAGenerated
            UTP_rRNA_cost += UTP_count*RNAGenerated
            CTP_rRNA_cost += CTP_count*RNAGenerated
            GTP_rRNA_cost += GTP_count*RNAGenerated
            ATP_trsc_energy_cost += (ATP_count + UTP_count + GTP_count + CTP_count)*RNAGenerated

        elif locusDic['Type'] == 'tRNA':
            Produced_RNA = 'Produced_R_' +locusNum
            RNAGenerated = int(countsDic[Produced_RNA][-1] - countsDic[Produced_RNA][-2])


            rnasequence = locusDic['RNAsequence']
            ATP_count = rnasequence.count('A')
            UTP_count = rnasequence.count('U')
            CTP_count = rnasequence.count('C')
            GTP_count = rnasequence.count('G')

            ATP_tRNA_cost += ATP_count*RNAGenerated
            UTP_tRNA_cost += UTP_count*RNAGenerated
            CTP_tRNA_cost += CTP_count*RNAGenerated
            GTP_tRNA_cost += GTP_count*RNAGenerated
            ATP_trsc_energy_cost += (ATP_count + UTP_count + GTP_count + CTP_count)*RNAGenerated

    countsDic['ATP_mRNA_cost'].append(ATP_mRNA_cost);countsDic['ATP_rRNA_cost'].append(ATP_rRNA_cost);countsDic['ATP_tRNA_cost'].append(ATP_tRNA_cost)
    countsDic['UTP_mRNA_cost'].append(UTP_mRNA_cost); countsDic['UTP_rRNA_cost'].append(UTP_rRNA_cost); countsDic['UTP_tRNA_cost'].append(UTP_tRNA_cost)
    countsDic['GTP_mRNA_cost'].append(GTP_mRNA_cost); countsDic['GTP_rRNA_cost'].append(GTP_rRNA_cost); countsDic['GTP_tRNA_cost'].append(GTP_tRNA_cost)
    countsDic['CTP_mRNA_cost'].append(CTP_mRNA_cost); countsDic['CTP_rRNA_cost'].append(CTP_rRNA_cost); countsDic['CTP_tRNA_cost'].append(CTP_tRNA_cost)

    countsDic['ATP_trsc_cost'].append(ATP_trsc_energy_cost)

    return None


def calculatetRNAchargingCosts(sim_properties):

    """
    
    
    Description: calculate the cost of ATP in CME tRNA charging 

    """

    countsDic = sim_properties['counts']
    
    ATP_tRNACharging_cost = countsDic['M_atp_c'][-2] - countsDic['M_atp_c'][-1]

    countsDic['ATP_tRNAcharging_cost'].append(ATP_tRNACharging_cost)

    # tRNAmap = sim_properties['trna_map']

    # for tRNAaa, tRNAlist in tRNAmap.items():
    #     # For single amino acid, multiple tRNAs can act as carriers.

    #     aaCostID = tRNAaa + '_cost'
        
    #     tRNANum = len(tRNAlist)
        
    #     perchargedtRNACost = int(countsDic[aaCostID][-1]/tRNANum)

    #     for tRNA in tRNAlist:
    #         chargedtRNAID = tRNA + '_ch'
    #         if countsDic[chargedtRNAID][-1] >= perchargedtRNACost:

    #             countsDic[chargedtRNAID][-1] = countsDic[chargedtRNAID][-1] - perchargedtRNACost

    #         else:
    #             # print(chargedtRNAID + ' runs out at time ' + str(sim_properties['time_second'][-1]))
    #             countsDic[chargedtRNAID][-1] = 0


    return None


def calculateDegradationCosts(sim_properties):

    """
    Input: sim_properties Dictionary

    Return: None

    Called by hookSimulation

    Description: Calculate the cost of ATP and productions of NMPs in degradation of mRNAs and append costs into corresponding species
    
    1 ATP hydrolysis per bp to degradate the mRNA 
    
    """

    countsDic = sim_properties['counts']
    genome = sim_properties['genome']
    NMP_recycle_counters = {'AMP_mRNAdeg_cost':'M_amp_c', 'UMP_mRNAdeg_cost':'M_ump_c', 'CMP_mRNAdeg_cost':'M_cmp_c', 'GMP_mRNAdeg_cost':'M_gmp_c'}

    ATP_mRNAdeg_cost = 0

    AMP_mRNAdeg_cost = 0
    UMP_mRNAdeg_cost = 0
    CMP_mRNAdeg_cost = 0
    GMP_mRNAdeg_cost = 0

    for locusTag, locusDict in genome.items():
        
        if locusDict["Type"] == 'protein':
            locusNum = locusTag.split('_')[1]
            Degradated_mRNA = 'Degradated_mRNA_' + locusNum
            mRNADegradated = int(countsDic[Degradated_mRNA][-1]- countsDic[Degradated_mRNA][-2])

            rnasequence = locusDict['RNAsequence']
            AMP_count = rnasequence.count('A')
            UMP_count = rnasequence.count('U')
            CMP_count = rnasequence.count('C')
            GMP_count = rnasequence.count('G')


            AMP_mRNAdeg_cost += AMP_count*mRNADegradated
            UMP_mRNAdeg_cost += UMP_count*mRNADegradated
            CMP_mRNAdeg_cost += CMP_count*mRNADegradated
            GMP_mRNAdeg_cost += GMP_count*mRNADegradated

            ATP_mRNAdeg_cost += (AMP_count + UMP_count + CMP_count + GMP_count)*mRNADegradated

    
    countsDic['AMP_mRNAdeg_cost'].append(AMP_mRNAdeg_cost)
    countsDic['UMP_mRNAdeg_cost'].append(UMP_mRNAdeg_cost)
    countsDic['CMP_mRNAdeg_cost'].append(CMP_mRNAdeg_cost)
    countsDic['GMP_mRNAdeg_cost'].append(GMP_mRNAdeg_cost)

    
    countsDic['ATP_mRNAdeg_cost'].append(ATP_mRNAdeg_cost)



    return None



def calculateTranslationCosts(sim_properties):
    """
    Input: sim_properties dictionary

    Return: None

    Called by hookSimulation

    Description: Calculate the GTP and amino acids with charged tRNA costs of the tranlation reactions based on the number of Produced_P_XXXX 
    
    """
    countsDic = sim_properties['counts']

    genome = sim_properties['genome']
    GTP_translate_cost = 0

    aaCostMap = sim_properties['aaCostMap']
    

    aaCostCounts =   {'ALA_cost': 0,'ARG_cost': 0,'ASN_cost': 0,'ASP_cost': 0,'CYS_cost': 0,'GLU_cost': 0,'GLN_cost': 0,'GLY_cost': 0,'HIS_cost': 0,
    'ILE_cost': 0,'LEU_cost': 0,'LYS_cost': 0,'MET_cost': 0,'PHE_cost': 0,'PRO_cost': 0,'SER_cost': 0,'THR_cost': 0,'TRP_cost': 0,'TYR_cost': 0,
    'VAL_cost': 0}


    for locusTag, locusDict in genome.items():
        
        if locusDict["Type"] == 'protein':

            locusNum = locusTag.split('_')[1]
            Produced_Ptn = "Produced_P_" + locusNum
            proteinGenerated = int(countsDic[Produced_Ptn][-1] - countsDic[Produced_Ptn][-2])

            aasequence = locusDict["AAsequence"]

            GTP_translate_cost += proteinGenerated*len(aasequence)*2
            
            for aa, aaCostStr in aaCostMap.items():
                # Pat attention to do accumulation of keys in dictionary
                aaCostCounts[aaCostStr] += proteinGenerated*aasequence.count(aa)

    # for aaCostStr, aaCost in aaCostCounts.items():
    #     # print(aaCostStr, aaCost)
    #     countsDic[aaCostStr].append(aaCost)

    countsDic['GTP_translate_cost'].append(GTP_translate_cost)

    return None

def calculateTranslocationCosts(sim_properties):
    """
    Input: sim_properties dictionary

    Return: None

    Called by hookSimulation

    Description: calculate and append the ATP cost for translocation of membrane proteins 

    1 ATP per 10 AAs
    """


    memPtnsList = sim_properties['memPtnsList']

    countsDic = sim_properties['counts']

    genome = sim_properties['genome']
    
    translocateCost = 0

    for locusTag in memPtnsList:
        
        locusNum = locusTag.split('_')[1]
        
        ptnID = 'P_' + locusNum
        
        aasequence = genome[locusTag]['AAsequence']

        translocateCost +=  int((len(aasequence)/10)*(countsDic[ptnID][-1] - countsDic[ptnID][-2]))
                
    countsDic['ATP_transloc_cost'].append(translocateCost)

    return None
    

def communicateCostsToMetabolism(sim_properties):

    """
    Input: sim_properties dictionary

    Returns: None
    
    Called by: hookSimulation

    Description: Pass the recorded consumpution of nucleotides in different processes into the counts of metabolites (such as ATP_trsc_cost, ATP_mRNA_cost into M_atp_c)

    """

    # Species M_ppi_c, M_pi_c, M_atp_c, M_ctp_c, M_utp_c, M_gtp_c, M_adp_c, M_cdp_c, M_udp_c, M_gdp_c, M_datp_c, M_dttp_c, M_dctp_c, M_dgtp_c,
    # M_amp_c, M_ump_c, M_cmp_c, M_gmp_c 's tragectories are appended by one new element. Also true for nucleotide_accumulative.
    
    countsDic = sim_properties['counts']



    ###################################
    #####  NTP Costs ####
    ###################################
    # the cost of atp in tRNA charging is already represented in the CME reactions
    NTPsCostMap = {'M_atp_c':{'syn_cost':['ATP_mRNA_cost', 'ATP_tRNA_cost', 'ATP_rRNA_cost'], 'energy_cost': ['ATP_trsc_cost', 'ATP_mRNAdeg_cost', 'ATP_DNArep_cost']},
                    'M_gtp_c':{'syn_cost': ['GTP_mRNA_cost','GTP_tRNA_cost','GTP_rRNA_cost'], 'energy_cost':['GTP_translate_cost']},
                    'M_ctp_c':{'syn_cost': ['CTP_mRNA_cost','CTP_tRNA_cost','CTP_rRNA_cost']},
                    'M_utp_c':{'syn_cost': ['UTP_mRNA_cost','UTP_tRNA_cost','UTP_rRNA_cost']},
                    }

   
    # ppi is generated in synthesis reactions.
    ppiCost = 0

    # pi is generated in energy related reactions.
    piCost = 0

    for NTPID, subDic in NTPsCostMap.items():
        
        nucName = NTPID[-5]

        NDPID = 'M_' + nucName + 'dp_c'

        NTPAccumulativeID = NTPID + '_accumulative'
        
        NTPAccumulativeCount = countsDic[NTPAccumulativeID][-1]
        
        # Pure NTP cost within 1 second
        NTPCost_energy = 0

        NTPCost_syn = 0

        # Counts of NTP from last second
        NTPCounts = countsDic[NTPID][-1]
        
        for function, processes in subDic.items():
            # print(function, processes)
            
            if function == 'energy_cost':
                for process in processes:
                    
                    NTPCost_energy += countsDic[process][-1]
                    piCost += countsDic[process][-1]

            else:
                for process in processes:

                    NTPCost_syn += countsDic[process][-1]                    
                    ppiCost += countsDic[process][-1]

        NewNTPAccumulative = NTPCounts + NTPAccumulativeCount - NTPCost_energy - NTPCost_syn


        if NTPID == 'M_atp_c':
            # Since M_atp_c has been appended in updataCMEcounts functions
            countsDic[NTPID][-1] = max(NewNTPAccumulative,1)
            countsDic[NTPAccumulativeID].append(min(NewNTPAccumulative,0))
        else:
            countsDic[NTPID].append(max(NewNTPAccumulative,1))
            countsDic[NTPAccumulativeID].append(min(NewNTPAccumulative,0))

        # If the NTP shortage happens, how to partition the NTP cost in energetic reactions and synthetic reactions ? 
        # Weight based on their consumptions
            
        if NTPCost_energy + NTPCost_syn != 0:
            # NDPCount always nonnegative
            NDPCount = int(min(NTPCost_energy, NTPCost_energy + NTPCost_energy*NewNTPAccumulative/(NTPCost_energy + NTPCost_syn)))
            countsDic[NDPID].append(countsDic[NDPID][-1] + NDPCount)
            piCost += min(0, NTPCost_energy*NewNTPAccumulative/(NTPCost_energy + NTPCost_syn))
            ppiCost += min(0, NTPCost_syn*NewNTPAccumulative/(NTPCost_energy + NTPCost_syn))

        else:
            # No NTP cost so no NDP, pi, and ppi generated
            countsDic[NDPID].append(countsDic[NDPID][-1])
    ###################################
    #####  dNTP Costs ####
    ###################################


    dNTPsCostMap = {'M_datp_c':'dATP_DNArep_cost', 'M_dttp_c': 'dTTP_DNArep_cost','M_dctp_c':'dCTP_DNArep_cost', 'M_dgtp_c':'dGTP_DNArep_cost'}
    
    for dNTPID, process in dNTPsCostMap.items():
        
        dNTPAccumulativeID = dNTPID + '_accumulative'

        dNTPAccumulativeCount = countsDic[dNTPAccumulativeID][-1]

        dNTPCount = countsDic[dNTPID][-1]

        dNTPCost = countsDic[process][-1]

        newdNTPaccumulative = dNTPCount + dNTPAccumulativeCount - dNTPCost

        countsDic[dNTPID].append(max(newdNTPaccumulative, 1))

        countsDic[dNTPAccumulativeID].append(min(newdNTPaccumulative,0))

        ppiCost += dNTPCost + min(0, newdNTPaccumulative)
    
    ###################################
    #####  NMP Costs/generation ####
    ###################################

    AMPsCostMap = {'M_amp_c':'AMP_mRNAdeg_cost', 'M_ump_c':'UMP_mRNAdeg_cost', 'M_cmp_c': 'CMP_mRNAdeg_cost','M_gmp_c': 'GMP_mRNAdeg_cost'  }

    for AMPID, process in AMPsCostMap.items():

        AMPAccumulativeID = AMPID + '_accumulative'

        AMPAccumulativeCount = countsDic[AMPAccumulativeID][-1]

        AMPCount = countsDic[AMPID][-1]
        
        AMPCost = countsDic[process][-1]

        newAMPaccumulative = AMPCount + AMPAccumulativeCount + AMPCost
        
        if AMPID == 'M_amp_c':
            countsDic[AMPID][-1] = max(newAMPaccumulative, 1)
        else:
            countsDic[AMPID].append(max(newAMPaccumulative, 1))
        countsDic[AMPAccumulativeID].append(min(newAMPaccumulative, 0))
    
    countsDic['M_pi_c'].append(countsDic['M_pi_c'][-1] + piCost )
    countsDic['M_ppi_c'][-1] = countsDic['M_ppi_c'][-1] + ppiCost

        
    return None


#########################################################################################
def updateSA(sim_properties):
    """
    Input: sim_properties dictionary

    Return: None

    Description: Calculate and update the current cell membrane surface and volume
    """
    genome = sim_properties['genome']
    
    countsDic = sim_properties['counts']
    avgProtSA = 28 #24.75 # nm^2, average protein surface area to produce expected 47% coverage for 9.6K membrane proteins
    
#     lipidSizes = {
#     'M_clpn_c':0.4,
#     'M_chsterol_c':0.35, # 0.35, test for Chol. value smaller
#     'M_sm_c':0.45,
#     'M_pc_c':0.55,
#     'M_pg_c':0.6,
#     'M_galfur12dgr_c':0.6,
#     'M_fa_c':0.5, # Should this be here??
#     'M_12dgr_c':0.5, # Scale down, should cdp-dag be added???
#     'M_pa_c':0.5,
#     'M_cdpdag_c':0.5,
#     }
    
    lipidSizes = {
    'M_clpn_c':0.4,
    'M_chsterol_c':0.35, 
    'M_sm_c':0.45,
    'M_pc_c':0.55,
    'M_pg_c':0.6,
    'M_galfur12dgr_c':0.6,
    'M_12dgr_c':0.5, 
    'M_pa_c':0.5,
    'M_cdpdag_c':0.5,
    }
    # 0.4 nm^2
    
    lipidSA = 0
    
    for lipid, size in lipidSizes.items():
        
        lipidSA += countsDic[lipid][-1]*size
    
    # Lipid accounts for 0.513 weight of surface area.
    lipidSA = int(lipidSA*0.513)
    
    sim_properties['SA']['SA_lipid'].append(lipidSA)

    # 93 membrane proteins
    memPtnsList = ['JCVISYN3A_0005','JCVISYN3A_0008', 'JCVISYN3A_0009', 'JCVISYN3A_0010', 'JCVISYN3A_0011', 'JCVISYN3A_0030', 
                   'JCVISYN3A_0034', 'JCVISYN3A_0060','JCVISYN3A_0095', 'JCVISYN3A_0113','JCVISYN3A_0114','JCVISYN3A_0116','JCVISYN3A_0117',
                   'JCVISYN3A_0132', 'JCVISYN3A_0143','JCVISYN3A_0146','JCVISYN3A_0164','JCVISYN3A_0165', 'JCVISYN3A_0166', 'JCVISYN3A_0167', 
                   'JCVISYN3A_0168', 'JCVISYN3A_0169', 'JCVISYN3A_0195', 'JCVISYN3A_0196', 'JCVISYN3A_0197','JCVISYN3A_0235','JCVISYN3A_0239',
                   'JCVISYN3A_0248','JCVISYN3A_0249','JCVISYN3A_0296','JCVISYN3A_0304','JCVISYN3A_0314','JCVISYN3A_0317','JCVISYN3A_0326',
                   'JCVISYN3A_0332','JCVISYN3A_0338','JCVISYN3A_0345', 'JCVISYN3A_0346','JCVISYN3A_0371','JCVISYN3A_0372','JCVISYN3A_0379',
                   'JCVISYN3A_0388','JCVISYN3A_0398','JCVISYN3A_0399','JCVISYN3A_0411','JCVISYN3A_0425', 'JCVISYN3A_0426', 'JCVISYN3A_0427', 
                   'JCVISYN3A_0428','JCVISYN3A_0439','JCVISYN3A_0440','JCVISYN3A_0478','JCVISYN3A_0481','JCVISYN3A_0505','JCVISYN3A_0516',
                   'JCVISYN3A_0601','JCVISYN3A_0639', 'JCVISYN3A_0641', 'JCVISYN3A_0642', 'JCVISYN3A_0643', 'JCVISYN3A_0652', 'JCVISYN3A_0685', 
                   'JCVISYN3A_0686', 'JCVISYN3A_0691','JCVISYN3A_0696', 'JCVISYN3A_0706', 'JCVISYN3A_0707', 'JCVISYN3A_0708', 'JCVISYN3A_0774', 
                   'JCVISYN3A_0777','JCVISYN3A_0778','JCVISYN3A_0779', 'JCVISYN3A_0787', 'JCVISYN3A_0789', 'JCVISYN3A_0790', 'JCVISYN3A_0791', 
                   'JCVISYN3A_0792','JCVISYN3A_0795', 'JCVISYN3A_0797', 'JCVISYN3A_0822', 'JCVISYN3A_0827', 'JCVISYN3A_0830', 'JCVISYN3A_0835',
                   'JCVISYN3A_0836', 'JCVISYN3A_0839', 'JCVISYN3A_0852','JCVISYN3A_0870', 'JCVISYN3A_0872', 'JCVISYN3A_0876', 'JCVISYN3A_0878', 
                   'JCVISYN3A_0879', 'JCVISYN3A_0881','JCVISYN3A_0908']
    
    memPtnCnt = 0
    

    for locusTag in memPtnsList:
        
        locusNum = locusTag.split('_')[1]
        
        ptnID = 'P_' + locusNum
        
        memPtnCnt += countsDic[ptnID][-1]

#     memPtnCnt = memPtnCnt + sim_properties['counts']['Degradosome']
        
    ptnSA = int(memPtnCnt*avgProtSA)
    
    sim_properties['SA']['SA_ptn'].append(ptnSA)
    
    sim_properties['SA']['SA_nm2'].append(int(lipidSA + ptnSA))
    
    sim_properties['SA']['SA_m2'].append(int(lipidSA + ptnSA)/1e18) # Unit m^2
    
    radius_2V = 2.502e-7   # m

    cyto_radius_nm_equivalent_sphere = np.sqrt(sim_properties['SA']['SA_m2'][-1]/(4*np.pi))   #m
    
    cyto_radius_nm = min(cyto_radius_nm_equivalent_sphere, radius_2V)


    sim_properties['volume_L'].append((4/3)*np.pi*(cyto_radius_nm)**3*1000)
    
    if cyto_radius_nm_equivalent_sphere > radius_2V:
        
        sim_properties['division_started'].append(True)
    
    # print('SA: ', sim_properties['SA'])
    # print('V: ', sim_properties['volume'])
    
    # print('cyto radius: ', sim_properties['cyto_radius'])
    # print('cyto radius nm: ', sim_properties['cyto_radius_nm'])
    
    return None



def saveConc(sim_properties,odeResults, odemodel):
    """
    Description: save the concentration of metabolites in ODE into sim_properties['conc']
    """     

    concDict = sim_properties['conc']

    resFinal = odeResults[-1,:]
    
    metIDDict = odemodel.getMetDict()

    for metID, Num in metIDDict.items():
        concDict[metID].append(resFinal[Num])

    return None

