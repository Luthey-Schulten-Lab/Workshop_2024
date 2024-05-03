"""
Author: Enguang Fu

Date: March 2024

Based on Zane Thornburg's Script

General functions that create genetic information processing reactions for individual genes.
Only define the species and reactions
Particle numbers are given in initializeCounts function

"""
# import numpy as np
import pandas as pd
import replication
import GIP_rates


#########################################################################################

def addGeneticInformationReactions(sim, sim_properties, gbfile):
    """
    Input: sim, sim_properties

    Return: None

    Called when restart new CME simulation

    Description: Define CME GIP species and add reactions; Create the list of species their initial counts of different GIP processes
    """
    dna = gbfile
    
    genome = sim_properties['genome']

    # add chromosome replication initiation and replication reactions
    # Genes are defined in addReplication reaction.
    # Numbers of species involved in initiation are already given in this function;
    # For replication, the numbers of intergenic species are defined, numbers of genes are not.

    ini_list, ini_counts = replication.addRepInit(sim, sim_properties)

    # 493 genes replicated G_XXXX, Produced_G_XXXX, JCVISYN3A_XXXX_inter
    rep_list, rep_counts = replication.addReplication(sim,sim_properties, gbfile = dna)
    
    # print('Replication Reactions Initialized')

    # define secy protein that is used translocation reactions for membrane proteins
    sim.defineSpecies(['P_0652','C_P_0652', 'Ribosome_mRNA_0652','Produced_P_0652'])
    

    currenttime_second = sim_properties['time_second'][-1]


    if currenttime_second == 0:
        flagListini = True
        
        proteomics = pd.read_excel('../input_data/initial_concentrations.xlsx', sheet_name='Comparative Proteomics')

        # 493 gene all have transcription reactions RNAP_G_XXXX,  R_XXXX, and Produced_R_XXXX
        trsc_list = []; trsc_counts = []
        # 455 gene coding for proteins: P_XXXX, Ribosome_mRNA_XXXX, Produced_P_XXXX
        translation_list = []; translation_counts = []
        # 455 gene coding for proteins: Degradosome_mRNA_XXXX, Degradated_mRNA_XXXX
        Deg_list = []; Deg_counts = []
    else:
        flagListini = False

    # define LEU_cost, LEU_cost_unpaid in CME simulation beforehand, since translation and tRNA charging will use them as reactants
    unpaidaaCostMap = sim_properties['unpaidaaCostMap']
    aaCostMap = sim_properties['aaCostMap']

    for unpaidaaCost in unpaidaaCostMap.values():
        sim.defineSpecies([unpaidaaCost])
        
        
    for aaCost in aaCostMap.values():
        sim.defineSpecies([aaCost])
    

    memPtnsList = sim_properties['memPtnsList']

# 'Type' = ['protein', 'ncRNA', 'gene', 'rRNA', 'tRNA', 'tmRNA']

    for locusTag, locusDict in genome.items():
        locusNum = locusTag.split('_')[1]
        rnasequence = locusDict["RNAsequence"]

        if locusDict['Type'] == 'gene':
            # Nothing to three pseudo genes 
            # print(locusTag)
            None

        
        elif locusDict['Type'] == 'rRNA':
            if (len(rnasequence) > 2*sim_properties['rnap_spacing']):
            # Gene_0068, Gene_0069, Gene_0533 and Gene_0534
                transcriptionLong(sim, sim_properties, locusTag, rnasequence)
                if flagListini:
                    trsc_XXXX = ['RNAP_G_'+locusNum, 'R_'+ locusNum, 'Produced_R_'+locusNum]
                    
                    trsc_list.extend(trsc_XXXX)
                    
                    trsc_counts.extend([0,1,0])

            else: 
                transcription(sim, sim_properties, locusTag, rnasequence)
                if flagListini:
                    trsc_XXXX = ['RNAP_G_'+locusNum, 'R_'+ locusNum, 'Produced_R_'+locusNum]
                    
                    trsc_list.extend(trsc_XXXX)
                    
                    trsc_counts.extend([0,1,0])

        elif locusDict["Type"] == 'protein':
            # Add ribosome bindind and translation of each mRNA
            transcription(sim, sim_properties, locusTag, rnasequence)
            # Add degradosome binding and degradatation of each mRNA
            degradation_mrna(sim, sim_properties, locusNum, rnasequence)
            if flagListini:
                trsc_XXXX = ['RNAP_G_'+locusNum, 'R_'+ locusNum, 'Produced_R_'+locusNum]
                
                trsc_list.extend(trsc_XXXX)
                
                trsc_counts.extend([0,1,0]) 

                deg_XXXX = ['Degradosome_mRNA_' + locusNum, 'Degradated_mRNA_' + locusNum]
                Deg_list.extend(deg_XXXX)
                Deg_counts.extend([0,0])

            aasequence = locusDict["AAsequence"]

            if locusTag in memPtnsList:
                
                membranePtnTranslationWithCost(sim, sim_properties, locusNum, aasequence)
                if flagListini:
                    tranlation_XXXX = [ 'P_' + locusNum,'C_P_' + locusNum, 'Ribosome_mRNA_'+locusNum, 'Produced_P_' + locusNum]
                    PtnIniCounts = proteomics.loc[proteomics["Locus Tag"] == locusTag]['Sim. Initial Ptn Cnt'].values[0]
                    translation_list.extend(tranlation_XXXX)
                    translation_counts.extend([PtnIniCounts, 0, 0, 0])
                    
            else:
            
                translationWithCost(sim, sim_properties, locusNum, aasequence)
                if flagListini:
                    tranlation_XXXX = [ 'P_' + locusNum, 'Ribosome_mRNA_'+locusNum, 'Produced_P_' + locusNum]
                    PtnIniCounts = proteomics.loc[proteomics["Locus Tag"] == locusTag]['Sim. Initial Ptn Cnt'].values[0]
        
                    translation_list.extend(tranlation_XXXX)
                    translation_counts.extend([PtnIniCounts, 0, 0])

        elif locusDict["Type"] == 'tRNA':
            # 29 tRNAs
            transcription(sim, sim_properties, locusTag, rnasequence) 
            # R_XXXX means free tRNA; acoording to cell 2022, we assume 40 free tRNA, 160 charged tRNA at initial conditions
            if flagListini:
                trsc_XXXX = ['RNAP_G_'+locusNum, 'R_'+ locusNum, 'Produced_R_'+locusNum]
                
                trsc_list.extend(trsc_XXXX)
                
                trsc_counts.extend([0,40,0])
            
        else:
            # print(locusTag)
            # tmRNA G_0158 and ncRNA G_0049 and G_0356
            transcription(sim, sim_properties, locusTag, rnasequence)
            if flagListini:
                trsc_XXXX = ['RNAP_G_'+locusNum, 'R_'+ locusNum, 'Produced_R_'+locusNum]
                
                trsc_list.extend(trsc_XXXX)
                
                trsc_counts.extend([0,1,0]) 

    tRNA_list, tRNA_counts = tRNAcharging(sim,sim_properties)  

    sim.defineSpecies(['RNAP_released'])
    sim.defineSpecies(['ribosome_released'])  
    sim.defineSpecies(['Degradosome_released'])


    if flagListini:
        trsc_list.append('RNAP_released');trsc_counts.append(0)
        translation_list.append('ribosome_released'); translation_counts.append(0)
        Deg_list.append('Degradosome_released'); Deg_counts.append(0)
        sim_properties['ini_list'] = ini_list; sim_properties['ini_counts'] = ini_counts
        sim_properties['rep_list'] = rep_list; sim_properties['rep_counts'] = rep_counts
        sim_properties['trsc_list'] = trsc_list; sim_properties['trsc_counts'] = trsc_counts
        sim_properties['translation_list'] = translation_list; sim_properties['translation_counts'] = translation_counts
        sim_properties['Deg_list'] = Deg_list; sim_properties['Deg_counts'] = Deg_counts
        sim_properties['tRNA_list'] = tRNA_list; sim_properties['tRNA_counts'] = tRNA_counts
        #print('The lists of GIP species and counts are passed to sim_properites Dictionary')
        
    #print('Transcription and translation reactions Defined')
    return None     

#########################################################################################

    

#########################################################################################
def transcription(sim, sim_properties, locusTag, rnasequence):

    locusNum = locusTag.split('_')[1]
    
    gene = 'G_' + locusNum

    rnaID  = 'R_' + locusNum
    RNAP_gene = 'RNAP_G_' + locusNum   #
    RNAProduced = 'Produced_R_' + locusNum   # The accumulative number of the newly synthesized mRNA from transcription
    RNAPReleased = 'RNAP_released' 
    species = [RNAP_gene, rnaID,  RNAProduced]
    sim.defineSpecies(species)

    # RNAP_binding is the function to calculate the rate of binding
    # RNAP_binding is independent with metabolism
    ########################################################################
    # Calculate RNAP binding rate that is dependent on promoter strength and cellvolume
    proxyPromoterStrength = sim_properties['promoters'][locusTag]
    # The volume of E coli is constant
    Ecoli_V = 1e-15  
    avgdr = 6.022e23

    cellVolume = sim_properties['volume_L'][-1]

    #binding_rate = 4*(180/765)*(proxyPromoterStrength/180)*Ecoli_V*avgdr/11400/60
    RNAP_binding_rate = 4*(180/765)*(proxyPromoterStrength/180)*Ecoli_V*avgdr/11400/60/avgdr/cellVolume

    ###############################################################################

    RNApol = 'RNAP'
    sim.addReaction((gene, RNApol), RNAP_gene, RNAP_binding_rate) 

    #RNAP_gene = 'RP_' + locusNum + '_C1'  
    ## new_RNA_ID = 'RP_' + locusNum + '_f' + '_C1' 
    ### Already CME ##
    # Transcription rate is dependent with the counts/concentrations of NTPs in the metabolism
    # The concentrations of NTPs are stored under sim_properties['counts]

    sim.addReaction(RNAP_gene, tuple([RNApol,gene,rnaID, RNAProduced, RNAPReleased]),GIP_rates.TranscriptionRate(sim_properties, locusTag, rnasequence))

    return None
#########################################################################################

#########################################################################################
def transcriptionLong(sim, sim_properties, locusTag, rnasequence):

    locusNum = locusTag.split('_')[1]
    
    gene = 'G_' + locusNum

    rnaID  = 'R_' + locusNum
    RNAP_gene = 'RNAP_G_' + locusNum   #
    RNAProduced = 'Produced_R_' + locusNum   # The accumulative number of the newly synthesized mRNA from transcription 
    RNAPReleased = 'RNAP_released' 
    species = [RNAP_gene,rnaID, RNAProduced]
    sim.defineSpecies(species)

    # RNAP_binding is the function to calculate the rate of binding
    # RNAP_binding is independent with metabolism
    ########################################################################
    # Calculate RNAP binding rate that is dependent on promoter strength and cellvolume
    proxyPromoterStrength = sim_properties['promoters'][locusTag]
    # The volume of E coli is constant
    Ecoli_V = 1e-15
    avgdr = 6.022e23

    cellVolume = sim_properties['volume_L'][-1]

    #binding_rate = 4*(180/765)*(proxyPromoterStrength/180)*Ecoli_V*avgdr/11400/60
    RNAP_binding_rate = 4*(180/765)*(proxyPromoterStrength/180)*Ecoli_V*avgdr/11400/60/avgdr/cellVolume

    ###############################################################################

    RNApol = 'RNAP'
    sim.addReaction((gene, RNApol), RNAP_gene, RNAP_binding_rate) 

    #RNAP_gene = 'RP_' + locusNum + '_C1'  
    ## new_RNA_ID = 'RP_' + locusNum + '_f' + '_C1' 
    ### Already CME ##
    # Transcription rate is dependent with the counts/concentrations of NTPs in the metabolism
    # The concentrations of NTPs are stored under sim_properties['counts]

    sim.addReaction(RNAP_gene, tuple([RNApol,gene,rnaID, RNAProduced, RNAPReleased]),GIP_rates.TranscriptionRate(sim_properties, locusTag, rnasequence))


    return None


#########################################################################################


#########################################################################################
def degradation_mrna(sim, sim_properties, locusNum, rnasequence): 

    # Use a two-step explicit binding model 
    mRNA = 'R_' + locusNum
    degradosome = 'Degradosome'
    # Deg_mRNA means Degradosome-mRNA
    Deg_mRNA = 'Degradosome_mRNA_' + locusNum
    Degradated_mRNA = 'Degradated_mRNA_' + locusNum
    species = [Deg_mRNA, Degradated_mRNA]
    Degradosome_released = 'Degradosome_released'
    sim.defineSpecies(species)
    
    cellVolume = sim_properties['volume_L'][-1]

    Ecoli_V = 1e-15
    avgdr = 6.022e23
    deg_bind_rate = 11*avgdr*Ecoli_V/60/7800/avgdr/cellVolume #2400 #7800 #1/M/s 0.00228
    # 11 # of cleavage events of RNase E in E.coli per minute per RNase E; 
    # 2400 and 7800 are the number of mRNA in slow and fast growing E.coli 
    # deg_bind_rate is the event that mRNA binds with degradosome per RNase E per mRNA per second in cellVolume
    # rate of binding with degradosome and mRNA is constant for all mRNAs

    sim.addReaction((mRNA, degradosome),Deg_mRNA, deg_bind_rate)
    
   
    
    # sim_properties['counts']['DM_' + locusNum] = 0  for communications between methods
    # mRNA degradation rate depends only on its length
    rna_deg_rate = GIP_rates.mrnaDegradationRate(rnasequence)
    
    sim.addReaction(Deg_mRNA, tuple([degradosome, Degradated_mRNA, Degradosome_released]), rna_deg_rate)

    return None
#########################################################################################


#########################################################################################
def translationWithCost(sim, sim_properties, locusNum, aasequence):
    """
    
    Description: translation reactions without translocation; unpaid aa cost defined on the product side
    """

    ptnID = 'P_' + locusNum
    mRNA = 'R_'+locusNum
    ribo_part = 'ribosomeP'
    mRNA_ribo = 'Ribosome_mRNA_'+locusNum
    PtnProduced = 'Produced_P_' + locusNum
    ribo_released = 'ribosome_released'
   
    species = [ptnID, mRNA_ribo, PtnProduced]
    sim.defineSpecies(species)

    cellVolume = sim_properties['volume_L'][-1]

    Ecoli_V = 1e-15
    avgdr = 6.022e23

    ribo_bind = 40*Ecoli_V*avgdr/60/6800/avgdr/cellVolume  # 0.002925642336248077 1/s
    # ribo_bind is the event that mRNA binds with ribosome per ribosome per mRNA per second in cellVolume
    # rate of binding between ribosome and mRNA is constant for all mRNAs

    sim.addReaction((mRNA, ribo_part), mRNA_ribo, ribo_bind)

   
    translation_rate = GIP_rates.TranslationRate(sim_properties, aasequence)
    
    translationProduct = [mRNA,ribo_part,ptnID, PtnProduced,ribo_released]
    unpaidaaCostMap = sim_properties['unpaidaaCostMap']

    # append LEU_cost_unpaid to the product side
    for aa in aasequence:
        if aa != '*':
            translationProduct.append(unpaidaaCostMap[aa])

    sim.addReaction(mRNA_ribo, tuple(translationProduct), translation_rate)

    return None
#########################################################################################


#########################################################################################

def membranePtnTranslationWithCost(sim, sim_properties, locusNum, aasequence):
    """
    
    Description: translation reactions with translocation; unpaid aa cost defined on the product side

    """


    ptnID = 'P_' + locusNum
    cyto_ptnID = 'C_P_' + locusNum
    mRNA = 'R_'+locusNum
    ribo_part = 'ribosomeP'
    mRNA_ribo = 'Ribosome_mRNA_'+locusNum
    PtnProduced = 'Produced_P_' + locusNum
    ribo_released = 'ribosome_released'

     # Avoid define P_0652 twice
    if locusNum != '0652':
        species = [ptnID, cyto_ptnID, mRNA_ribo, PtnProduced]
        sim.defineSpecies(species)

    cellVolume = sim_properties['volume_L'][-1]

    Ecoli_V = 1e-15
    avgdr = 6.022e23

    ribo_bind = 40*Ecoli_V*avgdr/60/6800/avgdr/cellVolume  # 0.002925642336248077 1/s
    # ribo_bind is the event that mRNA binds with ribosome per ribosome per mRNA per second in cellVolume
    # rate of binding between ribosome and mRNA is constant for all mRNAs

    sim.addReaction((mRNA, ribo_part), mRNA_ribo, ribo_bind)

    translation_rate = GIP_rates.TranslationRate(sim_properties, aasequence)

    translationProduct = [mRNA,ribo_part,ptnID, PtnProduced,ribo_released]
    unpaidaaCostMap = sim_properties['unpaidaaCostMap']

    # append LEU_cost_unpaid to the product side
    for aa in aasequence:
        if aa != '*':
            translationProduct.append(unpaidaaCostMap[aa])


    sim.addReaction(mRNA_ribo, (mRNA, ribo_part, cyto_ptnID, PtnProduced), translation_rate)




    translocation_rate = GIP_rates.TranslocationRate(aasequence)

    # Secy insert membrane proteins into membranes   
    secy = 'P_0652'

    sim.addReaction((cyto_ptnID, secy), (ptnID, secy), translocation_rate)

    return None

#########################################################################################


#########################################################################################
def translation(sim, sim_properties, locusNum, aasequence):
    """
    
    Description: translation reactions without translocation; no unpaid aa cost defined on the product side
    """

# Use a two-step explicit binding model 
    
    ptnID = 'P_' + locusNum
    mRNA = 'R_'+locusNum
    ribo_part = 'ribosomeP'
    mRNA_ribo = 'Ribosome_mRNA_'+locusNum
    PtnProduced = 'Produced_P_' + locusNum
    ribo_released = 'ribosome_released'
   
    species = [ptnID, mRNA_ribo, PtnProduced]
    sim.defineSpecies(species)

    cellVolume = sim_properties['volume_L'][-1]

    Ecoli_V = 1e-15
    avgdr = 6.022e23

    ribo_bind = 40*Ecoli_V*avgdr/60/6800/avgdr/cellVolume  # 0.002925642336248077 1/s
    # ribo_bind is the event that mRNA binds with ribosome per ribosome per mRNA per second in cellVolume
    # rate of binding between ribosome and mRNA is constant for all mRNAs

    sim.addReaction((mRNA, ribo_part), mRNA_ribo, ribo_bind)

   
    translation_rate = GIP_rates.TranslationRate(sim_properties, aasequence)
    
    sim.addReaction(mRNA_ribo, (mRNA,ribo_part,ptnID, PtnProduced,ribo_released), translation_rate)
    
    return None




#########################################################################################

def membranePtnTranslation(sim, sim_properties, locusNum, aasequence):


    ptnID = 'P_' + locusNum
    cyto_ptnID = 'C_P_' + locusNum
    mRNA = 'R_'+locusNum
    ribo_part = 'ribosomeP'
    mRNA_ribo = 'Ribosome_mRNA_'+locusNum
    PtnProduced = 'Produced_P_' + locusNum
    
     # Avoid define P_0652 twice
    if locusNum != '0652':
        species = [ptnID, cyto_ptnID, mRNA_ribo, PtnProduced]
        sim.defineSpecies(species)

    cellVolume = sim_properties['volume_L'][-1]

    Ecoli_V = 1e-15
    avgdr = 6.022e23

    ribo_bind = 40*Ecoli_V*avgdr/60/6800/avgdr/cellVolume  # 0.002925642336248077 1/s
    # ribo_bind is the event that mRNA binds with ribosome per ribosome per mRNA per second in cellVolume
    # rate of binding between ribosome and mRNA is constant for all mRNAs

    sim.addReaction((mRNA, ribo_part), mRNA_ribo, ribo_bind)

    translation_rate = GIP_rates.TranslationRate(sim_properties, aasequence)
    sim.addReaction(mRNA_ribo, (mRNA, ribo_part, cyto_ptnID, PtnProduced), translation_rate)

    translocation_rate = GIP_rates.TranslocationRate(aasequence)

    # Secy insert membrane proteins into membranes   
    secy = 'P_0652'

    sim.addReaction((cyto_ptnID, secy), (ptnID, secy), translocation_rate)

    return None

#########################################################################################


#########################################################################################
def tRNAcharging(sim, sim_properties):
    """
    Input: sim, sim_properties

    Return: tRNA_list:list of species in tRNA_charging, tRNA themselves not included, tRNA_counts: initial numbers of species

    Called when restarting CME simulation; Put behind the initiationMetabolites

    Description: define the species (except free tRNA and proteins) in tRNA charging reactions in CME and add reactions;
    create the lists of defined species and their initial counts lists
    
    """ 
    countsDic = sim_properties['counts']
    
    tRNA_list = ['M_atp_c', 'M_amp_c', 'M_ppi_c']; tRNA_counts = [countsDic['M_atp_c'][-1], countsDic['M_amp_c'][-1],countsDic['M_ppi_c'][-1]]
    sim.defineSpecies(tuple(tRNA_list))
    
    # parameters for charging of 20 amino acids
    RXNS_params = pd.read_excel( '../input_data/kinetic_params.xlsx', sheet_name='tRNA Charging')

    # Mapping LEU to tRNAs
    # trna_map {'LEU': ['R_0070', 'R_0423', 'R_0506'],...}

    for tRNA_aa, rnaIDlist in sim_properties['trna_map'].items():

        tRNA_XXXX = []; tRNA_counts_XXXX = []
#         if tRNA_aa == 'GLN':
            
#             print('Sad')

#         else:
        # rxnID: LEUTRS
        rxnID = tRNA_aa + 'TRS'
        aaCost_unpaid = tRNA_aa + '_cost_unpaid'
        aaCost = tRNA_aa + '_cost'

        rxn_params = RXNS_params.loc[ RXNS_params["Reaction Name"] == rxnID ]

        # aaID: M_ala_L_c
        aaID = rxn_params.loc[ rxn_params["Parameter Type"] == 'amino acid' ]["Value"].values[0]

        aacount = countsDic[aaID][-1]
    
        # P_0163 .... already defined in translations
        synthetaseID = rxn_params.loc[ rxn_params["Parameter Type"] == 'synthetase' ]["Value"].values[0]
        # P_0163_atp
        synthetaseAtpID = synthetaseID + '_atp'
        # P_p163_atp_aa
        synthetaseAaID = synthetaseAtpID + '_aa'

        tRNA_XXXX.extend([aaID, synthetaseAtpID, synthetaseAaID])
        tRNA_counts_XXXX.extend([aacount,0,0])

        sim.defineSpecies(tuple([aaID, synthetaseAtpID, synthetaseAaID]))


        sim.addReaction(tuple([synthetaseID, 'M_atp_c']), synthetaseAtpID, rxn_params.loc[ rxn_params["Parameter Type"] == 'k_atp' ]["Value"].values[0])

    
        sim.addReaction(tuple([synthetaseAtpID, aaID]), synthetaseAaID, rxn_params.loc[ rxn_params["Parameter Type"] == 'k_aa' ]["Value"].values[0])

        for rnaID in rnaIDlist:

            synthetaseTrnaID = synthetaseAaID + '_' + rnaID

            sim.addReaction(tuple([synthetaseAaID, rnaID]), synthetaseTrnaID, rxn_params.loc[ rxn_params["Parameter Type"] == 'k_tRNA' ]["Value"].values[0])

            chargedTrnaID = rnaID + '_ch'

            # Marker to know the number of produced charged tRNA
            Produced_chargedTrna = 'Produced_' + chargedTrnaID
            
            sim.defineSpecies(tuple([synthetaseTrnaID, chargedTrnaID, Produced_chargedTrna]))

            tRNA_XXXX.extend([synthetaseTrnaID, chargedTrnaID, Produced_chargedTrna])
            # The initial number of 29 charged tRNA is 160

            tRNA_counts_XXXX.extend([0, 160, 0])

            sim.addReaction(synthetaseTrnaID, tuple(['M_amp_c', 'M_ppi_c', synthetaseID, chargedTrnaID, Produced_chargedTrna]), rxn_params.loc[ rxn_params["Parameter Type"] == 'k_cat' ]["Value"].values[0])

            sim.addReaction(tuple([aaCost_unpaid, chargedTrnaID]), tuple([rnaID, aaCost]), 1e5)

        tRNA_list.extend(tRNA_XXXX); tRNA_counts.extend(tRNA_counts_XXXX)

    # print('tRNA')
    
    return tRNA_list, tRNA_counts
#########################################################################################



    
    
    
