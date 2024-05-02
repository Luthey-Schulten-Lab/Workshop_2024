"""
Author: Enguang Fu

Date: March 2024

initialize constants, time traces of counts, membranes 
"""


import pandas as pd
import numpy as np

#########################################################################################
# Mapping gene bank file into Dict
def mapDNA(genome):
    DNAmap = {}
    chromosome_length = int((genome.features[0].location.end+1)/10)
    ori_ter_rotation_factor = int(chromosome_length/2)
    print(ori_ter_rotation_factor)   
    for feature in genome.features:
        strand = feature.strand     
        if strand == 1:    
            start = int(feature.location.start/10)        
            if start<=ori_ter_rotation_factor:
                start = start + ori_ter_rotation_factor           
            elif start>ori_ter_rotation_factor:               
                start = start - ori_ter_rotation_factor
            end = int(feature.location.end/10)
            if end<=ori_ter_rotation_factor:             
                end = end + ori_ter_rotation_factor
            elif end>ori_ter_rotation_factor:
                end = end - ori_ter_rotation_factor
        elif strand == -1:      
            start = int(feature.location.end/10)        
            if start<=ori_ter_rotation_factor:
                start = start + ori_ter_rotation_factor
            elif start>ori_ter_rotation_factor:
                start = start - ori_ter_rotation_factor
            end = int(feature.location.start/10)
            if end<=ori_ter_rotation_factor:            
                end = end + ori_ter_rotation_factor        
            elif end>ori_ter_rotation_factor:
                end = end - ori_ter_rotation_factor
        if feature.type == 'CDS':
            
            if('protein_id' in feature.qualifiers.keys()):
                # Excluding 3 pseudo genes
                # Totoal coding genes are 455
                locusTag = feature.qualifiers['locus_tag'][0]
                DNAmap[locusTag] = {}
                DNAmap[locusTag]['Type'] = 'protein'
                DNAmap[locusTag]['startIndex'] = [int(start)]
                DNAmap[locusTag]['originalStart'] = int(start)
                DNAmap[locusTag]['endIndex'] = [int(end)]           
                DNAmap[locusTag]['originalEnd'] = int(end)
                DNAmap[locusTag]['RNAsequence'] = str(feature.location.extract(genome.seq).transcribe())
                DNAmap[locusTag]['AAsequence'] = str(feature.location.extract(genome.seq).transcribe().translate(table=4))            
                DNAmap[locusTag]['GeneName'] = str(feature.qualifiers['product'][0])     
        elif feature.type == 'tRNA':
            locusTag = feature.qualifiers['locus_tag'][0]       
            DNAmap[locusTag] = {} 
            DNAmap[locusTag]['Type'] = 'tRNA'
            DNAmap[locusTag]['startIndex'] = [int(start)]          
            DNAmap[locusTag]['originalStart'] = int(start)
            DNAmap[locusTag]['endIndex'] = [int(end)]
            DNAmap[locusTag]['originalEnd'] = int(end) 
            DNAmap[locusTag]['RNAsequence'] = str(feature.location.extract(genome.seq).transcribe())
            DNAmap[locusTag]['GeneName'] = str(feature.qualifiers['product'][0])      
        elif feature.type == 'rRNA':           
            locusTag = feature.qualifiers['locus_tag'][0]           
            DNAmap[locusTag] = {}          
            DNAmap[locusTag]['Type'] = 'rRNA'
            DNAmap[locusTag]['startIndex'] = [int(start)]
            DNAmap[locusTag]['originalStart'] = int(start)
            DNAmap[locusTag]['endIndex'] = [int(end)]
            DNAmap[locusTag]['originalEnd'] = int(end)          
            DNAmap[locusTag]['RNAsequence'] = str(feature.location.extract(genome.seq).transcribe())
            DNAmap[locusTag]['GeneName'] = str(feature.qualifiers['product'][0])
            
        else:
            # 2 misc binding features, 2 ncRNA genes, 1 tm RNA gene, 
            try:
                locusTag = feature.qualifiers['locus_tag'][0]           
                DNAmap[locusTag] = {}          
                DNAmap[locusTag]['Type'] = feature.type
                DNAmap[locusTag]['startIndex'] = [int(start)]
                DNAmap[locusTag]['originalStart'] = int(start)
                DNAmap[locusTag]['endIndex'] = [int(end)]
                DNAmap[locusTag]['originalEnd'] = int(end)          
                DNAmap[locusTag]['RNAsequence'] = str(feature.location.extract(genome.seq).transcribe())
                DNAmap[locusTag]['GeneName'] = str(feature.qualifiers['product'][0])
            except:
                continue
            
    return DNAmap, chromosome_length
#########################################################################################

#########################################################################################
def initializeCME(sim, restartNum, sim_properties):
    """
    
    Decription: Set up the simulation time and hook and write intervals for CME simulation
    """
    
    # Initialize every restart CME simulation object
    restartInterval = sim_properties['restartInterval']
    hookInterval = sim_properties['hookInterval']
   
    sim.setWriteInterval(hookInterval)

    sim.setHookInterval(hookInterval)

    sim.setSimulationTime(restartInterval)
    
    return sim
#########################################################################################



#########################################################################################
def initializeConstants(sim_properties, genome3A):
    """
    
    Description: create constants including lists or sub dictionaries in sim_properties. 
        Sub dictionaries: tRNAMap, promoter strengths, and LocusNumtoIndex
        Lists: membrane proteins list, pseodogenes list

    """

    # Convert genebank file syn3A.gb into a multilayer dictionanry genome
    sim_properties['genome'], sim_properties['genome_length'] = mapDNA(genome3A)


    initializetRNAMap(sim_properties)
    initializeLocusNumtoIndex(sim_properties, genome3A)
    initializePromoterStrengths(sim_properties)


    sim_properties['memPtnsList'] = ['JCVISYN3A_0005','JCVISYN3A_0008', 'JCVISYN3A_0009', 'JCVISYN3A_0010', 'JCVISYN3A_0011', 'JCVISYN3A_0030', 
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
    
    sim_properties['aaCostMap'] = {"A":"ALA_cost", "R":"ARG_cost", 
        "N":"ASN_cost", "D":"ASP_cost", "C":"CYS_cost", "E":"GLU_cost", "Q":"GLN_cost", "G":"GLY_cost", 
            "H":"HIS_cost", "I":"ILE_cost", "L":"LEU_cost", "K":"LYS_cost", "M":"MET_cost", "F":"PHE_cost", 
        "P":"PRO_cost", "S":"SER_cost", "T":"THR_cost", "W":"TRP_cost", "Y":"TYR_cost", "V":"VAL_cost"}

    sim_properties['aaCost_list'] = list(sim_properties['aaCostMap'].values())

    # Mapping aa sequence to unpaid aacost species in CME
    sim_properties['unpaidaaCostMap'] = {'A': 'ALA_cost_unpaid', 'R': 'ARG_cost_unpaid', 'N': 'ASN_cost_unpaid', 'D': 'ASP_cost_unpaid', 'C': 'CYS_cost_unpaid',
      'E': 'GLU_cost_unpaid', 'Q': 'GLN_cost_unpaid', 'G': 'GLY_cost_unpaid', 'H': 'HIS_cost_unpaid', 'I': 'ILE_cost_unpaid', 
      'L': 'LEU_cost_unpaid', 'K': 'LYS_cost_unpaid', 'M': 'MET_cost_unpaid', 'F': 'PHE_cost_unpaid', 'P': 'PRO_cost_unpaid', 
      'S': 'SER_cost_unpaid','T': 'THR_cost_unpaid', 'W': 'TRP_cost_unpaid', 'Y': 'TYR_cost_unpaid', 'V': 'VAL_cost_unpaid'}


    sim_properties['unpaidaaCost_list'] = list(sim_properties['unpaidaaCostMap'].values())


    sim_properties['aa_list'] = ['M_ala__L_c','M_arg__L_c', 'M_asn__L_c', 'M_asp__L_c', 'M_cys__L_c',
        'M_glu__L_c', 'M_gly_c', 'M_his__L_c', 'M_ile__L_c', 'M_leu__L_c',
        'M_lys__L_c', 'M_met__L_c', 'M_phe__L_c', 'M_pro__L_c', 'M_ser__L_c',
        'M_thr__L_c', 'M_trp__L_c', 'M_tyr__L_c', 'M_val__L_c', 'M_gln__L_c']
    
    sim_properties['rnap_spacing'] = int(400)

    sim_properties['pseudoGenes'] = ['JCVISYN3A_0051', 'JCVISYN3A_0546','JCVISYN3A_0602']

    r_cell = 2.0*(10**-7) # m
    cellVolume_init = (4*np.pi/3)*r_cell**3*1000 # L
    sim_properties['volume_L'] = [cellVolume_init]

    # sim_properties['counts'] is a dictionary that record the number of all species per hookInterval
    sim_properties['counts'] = {}

    # sim_properties['conc'] is a dictionary that record the concentrations of metabolites after each ODE run
    sim_properties['conc'] = {}



    return None
#########################################################################################

#########################################################################################
def initializePromoterStrengths(sim_properties):
    """
    Input: sim_properties

    Return: None

    Called at the beginning of the whole simulation once

    Description:
    """

    # Promotoer Strength maintain in the whole cell cycle
    # Used to modify the transcription rates of mRNA, tRNA, and rRNA
    # Ribosomal proteins have higher transcription rate
    # kcat_mod: tRNA 85, rRNA 85*2.5 
    proxyPromoterStrengths = {}
    
    genome = sim_properties['genome']   

    proteomics = pd.read_excel('../input_data/initial_concentrations.xlsx', sheet_name='Comparative Proteomics')

    for locusTag, locusDict in genome.items():       
        if locusDict["Type"] == 'protein':
            PtnName = proteomics.loc[ proteomics['Locus Tag'] == locusTag ]['Gene Product'].values[0]
            if 'S ribosomal' in PtnName:
                # Error here??
                # For ribosomal proteins, the min is 500              
                proxyPromoterStrengths[locusTag] = min(765, 500+PtnCount)               
            else:               
                PtnCount = proteomics.loc[ proteomics['Locus Tag'] == locusTag ]['Sim. Initial Ptn Cnt'].values[0]
                proxyPromoterStrengths[locusTag] = min(765, max(45, PtnCount))              
        elif locusDict["Type"] == 'tRNA':

            #  proxyPromoterStrengths[locusTag] = 180
            # Different from cell paper
            proxyPromoterStrengths[locusTag] = 765           
        elif locusDict["Type"] == 'rRNA':
            # Different from cell paper
            proxyPromoterStrengths[locusTag] = 765*2.5 
            
        else:
            
            proxyPromoterStrengths[locusTag] = 10
    sim_properties['promoters'] = proxyPromoterStrengths

    return None
#########################################################################################



#########################################################################################
def initializeConcDictionary(sim_properties, odemodel):
    """"
    Called after the first time ODE run

    Description: initialize the trajectories of metabolites in ODE 
    """
    concDict = sim_properties['conc']
    metIDDict = odemodel.getMetDict()
    # {'M_ACP_c': 0, 'M_ACP_R_c': 1, 'M_apoACP_c': 2, ...}

    for species in metIDDict.keys():
        concDict[species] = []

    return None

#########################################################################################



#########################################################################################
def initializeLocusNumtoIndex(sim_properties, gbfile):
    """
    

    Called at the very beginning of the simulation once

    Description: Set up sub dictionary sim_properties['LocusNumtoIndex'] where the values are the genes' locusNums and keys are the start and end position of each gene
    """

    dna = gbfile
    gene_secondhalf = []
    gene_firsthalf = []
    LocusNumtoIndex = {}

    endposition = 543086
    # the index of the end position of 0910 is 543086
    # The length of the whole chromosome is len(dna.seq) is 543379
    position = 0

    gene_list = []
    for i in range(len(dna.features)):
        if ('product' in dna.features[i].qualifiers.keys()):
            #print(i) # This first statement works
            #print(dna.features[i].qualifiers['product'])
            if dna.features[i].qualifiers['product'][0]:# Figure out how to sort out for ribosomal operons?
                #print(dna.features[i].qualifiers['product'])
                gene_list.append(i)

    for gene in gene_list:
        locusTag = dna.features[gene].qualifiers['locus_tag'][0]
        locusNum = locusTag.split('_')[1] 
        start =  dna.features[gene].location.start.real
        end  = dna.features[gene].location.end.real
        if start < len(dna.seq)/2:
            gene_firsthalf.append(gene)
            if start == 0:
                LocusNumtoIndex[locusNum] = [position, end]
                position = end
            
            else:
                LocusNumtoIndex[locusNum] = [position, end]
                position = end
        else:
            gene_secondhalf.append(gene)

    position = endposition
    gene_secondhalf.reverse()

    for gene in gene_secondhalf:
        locusTag = dna.features[gene].qualifiers['locus_tag'][0]
        locusNum = locusTag.split('_')[1] 
        start =  dna.features[gene].location.start.real
        end  = dna.features[gene].location.end.real
        if end == endposition:

            LocusNumtoIndex[locusNum] = [start, position]
            position = start
            
        else:
            LocusNumtoIndex[locusNum] = [start, position]
            position = start

    sim_properties['LocusNumtoIndex'] = LocusNumtoIndex

    return None

#########################################################################################




#########################################################################################
def initializetRNAMap(sim_properties):
    """
    Input: sim_properties

    Return: None

    Called at the beginning of the simulation

    Description: sim_properties['trna_map']
    Set up the sub dictionary with amino acids and correspoding tRNAs {'LEU': ['R_0070', 'R_0423', 'R_0506'],...}
    """

    genome = sim_properties['genome']

    tRNA_map = {}

    for locusTag, locusDict in genome.items():
        
        if locusDict['Type'] == "tRNA":
            
            locusNum = locusTag.split('_')[1]
            
            rnaID = 'R_'+locusNum

            rnaName = locusDict['GeneName'].split('-')[1].upper()
            
            if rnaName not in tRNA_map:
                
                tRNA_map[rnaName] = [rnaID]
                
            else:
                
                tRNA_map[rnaName].append(rnaID)
                
    sim_properties['trna_map'] = tRNA_map

    return None



#########################################################################################



#########################################################################################
def initializeCosts(sim_properties):
    """
    
    Called once at the beginning of the simulation
    
    Description: Initialize the nucleotides cost species in counts dictionary

    """

    countsDic = sim_properties['counts']

    # nucleotides cost as energetic reactions NTPs -> NDPs + pi; For tRNA charging to AMP and ppi
    energy_cost_counters = {'GTP_translate_cost':'g', 'ATP_trsc_cost':'a', 'ATP_mRNAdeg_cost':'a', 
                            'ATP_DNArep_cost':'a', 'ATP_transloc_cost':'a','ATP_tRNAcharging_cost':'a'}

    # nucleotides cost in polymerization processes NTPs -> NMPs + ppi
    nuc_costs = {'ATP_mRNA_cost':'M_atp_c', 'CTP_mRNA_cost':'M_ctp_c', 'UTP_mRNA_cost':'M_utp_c', 'GTP_mRNA_cost':'M_gtp_c',
                    'ATP_tRNA_cost':'M_atp_c', 'CTP_tRNA_cost':'M_ctp_c', 'UTP_tRNA_cost':'M_utp_c', 'GTP_tRNA_cost':'M_gtp_c',
                    'ATP_rRNA_cost':'M_atp_c', 'CTP_rRNA_cost':'M_ctp_c', 'UTP_rRNA_cost':'M_utp_c', 'GTP_rRNA_cost':'M_gtp_c',
                    'dATP_DNArep_cost':'M_datp_c', 'dTTP_DNArep_cost':'M_dttp_c', 'dCTP_DNArep_cost':'M_dctp_c', 'dGTP_DNArep_cost':'M_dgtp_c'}
    #
    NMP_recycle_counters = {'AMP_mRNAdeg_cost':'M_amp_c', 'UMP_mRNAdeg_cost':'M_ump_c', 'CMP_mRNAdeg_cost':'M_cmp_c', 'GMP_mRNAdeg_cost':'M_gmp_c'}



    for value, key in energy_cost_counters.items():
        countsDic[value] = [int(0)]
        # countsDic[value+'_accumulative'] = [int(0)]
    
    for value, key in nuc_costs.items():
        countsDic[value] = [int(0)]
        # countsDic[value+'_accumulative'] = [int(0)]

    for value, key in NMP_recycle_counters.items():
        countsDic[value] = [int(0)]
        # countsDic[value+'_accumulative'] = [int(0)]
    
    monomerlist = ['M_atp_c', 'M_utp_c', 'M_ctp_c', 'M_gtp_c','M_datp_c','M_dttp_c','M_dctp_c','M_dgtp_c','M_amp_c','M_ump_c','M_cmp_c','M_gmp_c']
    
    for metID in monomerlist:
        # nucName = monomerName.split('_')[1]
        nucAccumulative = metID + '_accumulative'
        countsDic[nucAccumulative] = [int(0)]
    
    return None 
#########################################################################################



#########################################################################################
def initializeMediumConcs(sim_properties):
    
    sim_properties['medium'] = {}
    
    sim_medium = pd.read_excel('../input_data/initial_concentrations.xlsx', sheet_name='Simulation Medium')
    
    for row, nutrient in sim_medium.iterrows():
        
        metID = 'M_' + nutrient['Met ID']
        
        # The concentration of medium species is constant and put in another subDic.
        sim_properties['medium'][metID] = nutrient['Conc (mM)']
        

    return None
#########################################################################################


#########################################################################################
def initializeMetabolitesCounts(sim_properties):
    """
    Input: sim_properties

    Return: None

    Description: Add the initial counts of metabolites into counts dictionary
    
    """


    metabolite_ic = pd.read_excel('../input_data/initial_concentrations.xlsx', sheet_name='Intracellular Metabolites')
    countsDic = sim_properties['counts']

    for row, metabolite in metabolite_ic.iterrows():
        
        metID = 'M_' + metabolite['Met ID']
        
        # The first element of metabolite's trajectory is the initial counts at 0 second
        countsDic[metID] = [mMtoPart(metabolite['Init Conc (mM)'], sim_properties)]
    
    # print('Metabolism Initialized')
    return None
#########################################################################################


#########################################################################################
def initializeProteinMetabolitesCounts(sim_properties):
    """

    Called once at the very beginning of the simulation

    Description: Define and add the initial counts of forms of four proteins 'P_0233', 'P_0694', 'P_0234', 'P_0779' in phosphorelay; The phospholated form starts from 0
    """

    phosphorelayPtns = ['P_0233', 'P_0694', 'P_0234', 'P_0779']

    phosphorelayPtnsMap = {'P_0233':352, 'P_0694':289, 'P_0234':313, 'P_0779':830}

    # phosphorelayPtnsCounts = [352, 289, 313, 830]

    data_file =  '../input_data/initial_concentrations.xlsx'
    
    ptnMets = pd.read_excel(data_file, sheet_name='protein_metabolites')
    
    for index, row in ptnMets.iterrows():
        PtnID = row['Protein']
        if PtnID in phosphorelayPtns:

            metabolites = row['Metabolite IDs'].split(',')
        
            for metID in metabolites:
                if metID == metabolites[0]:

                    sim_properties['counts'][metID] = [phosphorelayPtnsMap[PtnID]]

                else:

                    sim_properties['counts'][metID] = [0]
        

    return None
#########################################################################################


#########################################################################################
def addGeneticInformationSpeciesCounts(sim,sim_properties):
    """
    Input: sim, sim_properties

    Return: None

    Called when restart new CME simulation

    Description: Initialize the trajectory of species in sim_properties dictionary; Add species counts to new CME simulation; 
    """

    # Called when a new CME simulation is created
    # GRP stands for genes, RNAs and proteins
    # Initialize numbers of each gene as 1, of each mRNA as 1 at the very beginning of the simulation
    # The number of proteins are determined based on experiment

    addinitiationCounts(sim, sim_properties)
    addReplicationCounts(sim, sim_properties)
    addTranscriptionCounts(sim, sim_properties)
    addTranslationCounts(sim, sim_properties)
    addDegradationCounts(sim, sim_properties)
    addtRNAChargingCounts(sim, sim_properties)
    addaaCostCounts(sim, sim_properties)
    
    return None
#########################################################################################


#########################################################################################
def addinitiationCounts(sim, sim_properties):
    time_second = sim_properties['time_second']
    countsDic = sim_properties['counts']
    ini_list = sim_properties['ini_list']
    ini_counts = sim_properties['ini_counts']

    if time_second[-1] == 0:
        for i in range(len(ini_list)):
            sim.addParticles(species = ini_list[i], count = ini_counts[i])
            countsDic[ini_list[i]] = [ini_counts[i]]

    else:
        for i in range(len(ini_list)):
            sim.addParticles(species = ini_list[i], count = countsDic[ini_list[i]][-1])

    return None
#########################################################################################



#########################################################################################

def addReplicationCounts(sim, sim_properties):
    time_second = sim_properties['time_second']
    countsDic = sim_properties['counts']
    rep_list = sim_properties['rep_list']
    rep_counts = sim_properties['rep_counts']


    if time_second[-1] == 0:
        for i in range(len(rep_list)):
            sim.addParticles(species = rep_list[i], count = rep_counts[i])
            countsDic[rep_list[i]] = [rep_counts[i]]

    else:
        for i in range(len(rep_list)):
            sim.addParticles(species = rep_list[i], count = countsDic[rep_list[i]][-1])
    return None


#########################################################################################
def addTranscriptionCounts(sim, sim_properties):

    time_second = sim_properties['time_second']
    countsDic = sim_properties['counts']

    trsc_list = sim_properties['trsc_list']
    trsc_counts = sim_properties['trsc_counts']


    if time_second[-1] == 0:

        for i in range(len(trsc_list)):
            # print(trsc_list[i])
            sim.addParticles(species = trsc_list[i], count = trsc_counts[i])
            countsDic[trsc_list[i]] = [trsc_counts[i]]
        
        # print("Transcription Initialized")

    else:
        for i in range(len(trsc_list)):
            sim.addParticles(species = trsc_list[i], count = countsDic[trsc_list[i]][-1])

    return None


#########################################################################################



#########################################################################################
def addTranslationCounts(sim, sim_properties):
    
    time_second = sim_properties['time_second']
    countsDic = sim_properties['counts']

    translation_list = sim_properties['translation_list']
    translation_counts = sim_properties['translation_counts']
    # When defining the number of each species, we need to give the a conditional statement 
    # to differ the very beginning of the whole cycle and the restart simulation
    if time_second[-1] == 0:

        for i in range(len(translation_list)):
            sim.addParticles(species = translation_list[i], count = translation_counts[i])
            countsDic[translation_list[i]] = [translation_counts[i]]
        # print('Translation Initialized')
    else:
        
        # Pass the last element of trajectory recorded in coundsDic to next newly restart CME simulation
        for i in range(len(translation_list)):
            sim.addParticles(species = translation_list[i], count = countsDic[translation_list[i]][-1])

    #print('Translation Initialized')

    return None

#########################################################################################



#########################################################################################
def addDegradationCounts(sim, sim_properties):
    time_second = sim_properties['time_second']
    countsDic = sim_properties['counts']
    Deg_list = sim_properties['Deg_list']
    Deg_counts = sim_properties['Deg_counts']

    if time_second[-1] == 0:
        for i in range(len(Deg_list)):
            sim.addParticles(species = Deg_list[i], count = Deg_counts[i])
            countsDic[Deg_list[i]] = [Deg_counts[i]]
    else:
        for i in range(len(Deg_list)):
            sim.addParticles(species = Deg_list[i], count = countsDic[Deg_list[i]][-1])


    return None

#########################################################################################

#########################################################################################
def addtRNAChargingCounts(sim, sim_properties):
    """
    
    
    Description: addParticles to tRNA charging related species, including atp, amp, ppi, R_XXXX_ch
    """
    
    time_second = sim_properties['time_second']
    countsDic = sim_properties['counts']
    tRNA_list = sim_properties['tRNA_list'] 
    tRNA_counts = sim_properties['tRNA_counts']
    
    if time_second[-1] == 0:
        for i in range(len(tRNA_list)):
            sim.addParticles(species = tRNA_list[i], count = tRNA_counts[i])
            countsDic[tRNA_list[i]] = [tRNA_counts[i]]
    else:
        for i in range(len(tRNA_list)):
            sim.addParticles(species = tRNA_list[i], count = countsDic[tRNA_list[i]][-1])
    return None

#########################################################################################

def addaaCostCounts(sim,sim_properties):
    """
    
    Description: addParticles to aa_cost, aa_cost_unpaid species to CME; initialize aa_cost, aa_cost_unpaid in countsDic
    """

    time_second = sim_properties['time_second']
    countsDic = sim_properties['counts']
    aaCost_list = sim_properties['aaCost_list']
    unpaidaaCost_list = sim_properties['unpaidaaCost_list']

    if time_second[-1] == 0:
        for i in range(len(aaCost_list)):
            sim.addParticles(species = aaCost_list[i], count = 0)
            countsDic[aaCost_list[i]] = [0]
        for i in range(len(unpaidaaCost_list)):
            sim.addParticles(species = unpaidaaCost_list[i], count = 0)
            countsDic[unpaidaaCost_list[i]] = [0]

    else:
        for i in range(len(aaCost_list)):
            sim.addParticles(species = aaCost_list[i], count = countsDic[aaCost_list[i]][-1])
   
        for i in range(len(unpaidaaCost_list)):
            sim.addParticles(species = unpaidaaCost_list[i], count = countsDic[unpaidaaCost_list[i]][-1])

    return None

#########################################################################################
def addEnzymesCounts(sim, sim_properties):
    # Initialize the numbers of three enzymes when new CME simulation is created
    countsDic = sim_properties['counts']
    species = ['RNAP','ribosomeP','Degradosome']
    sim.defineSpecies(species)

    if sim_properties['time_second'][-1] == 0:
        # Need to define species
        # The numbers for the very beginning of whole cell cycle
        sim.addParticles(species = 'RNAP', count = 187)
        sim.addParticles(species = 'ribosomeP', count = 500)
        sim.addParticles(species = 'Degradosome', count = 120)

        # For each species, a list is created with initial value
        countsDic['RNAP'] = [int(187)]
        countsDic['ribosomeP'] = [int(500)]
        countsDic['Degradosome'] = [int(120)]

    else:
        # Read the last element of the list and pass to next restart simulation
        sim.addParticles(species = 'RNAP',count = countsDic['RNAP'][-1] )
        sim.addParticles(species = 'ribosomeP',count = countsDic['ribosomeP'][-1] )
        sim.addParticles(species = 'Degradosome',count = countsDic['Degradosome'][-1] )


#########################################################################################




#########################################################################################
def mMtoPart(conc, sim_properties):
    """
    Convert particle counts to mM concentrations for the ODE Solver

    Parameters:
    particles (int): The number of particles for a given chemical species

    Returns:
    conc (float): The concentration of the chemical species in mM
    """

    ### Constants
    NA = 6.022e23 # Avogadro's
    
    count = max(1, int(round((conc/1000)*NA*sim_properties['volume_L'][-1])))

    return count
#########################################################################################




#########################################################################################
def initializeMembrane(sim_properties):
    """
    Input: sim_properties dictionary

    Return: None

    Description: Calculate the membrane surface based on the counts of lipids and membrane proteins P_XXXX and initialize the counts in sim_properties
    Calcualte the cellVolume based on 200 nm radius
    
    """

    
    # Calculate the surface area at the initial time

    genome = sim_properties['genome']
    
    countsDic = sim_properties['counts']

    avgProtSA = 28 #24.75 # nm^2, average protein surface area to produce expected 47% coverage for 9.6K membrane proteins
    
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
        
        lipidSA = lipidSA + countsDic[lipid][-1]*size
    
    # Lipid accounts for 0.513 weight of surface area.
    lipidSA = int(lipidSA*0.513)

    sim_properties['SA'] = {}

    sim_properties['SA']['SA_lipid'] = [lipidSA]

    
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
        
        memPtnCnt = memPtnCnt + countsDic[ptnID][-1]

    ptnSA = int(memPtnCnt*avgProtSA)

    sim_properties['SA']['SA_ptn'] = [ptnSA]
    
    sim_properties['SA']['SA_nm2'] = [int(lipidSA + ptnSA)]
    
    sim_properties['SA']['SA_m2'] = [int(lipidSA + ptnSA)/1e18]    # m^2

    radius_2V = 2.502e-7   # m

    cyto_radius_nm_equivalent_sphere = np.sqrt(sim_properties['SA']['SA_m2'][-1]/(4*np.pi))   #m

    sim_properties['division_started'] = [False]



    return None