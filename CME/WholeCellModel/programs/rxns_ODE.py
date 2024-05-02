"""
Authors: Zane Thornburg

A file to define all of the reactions in the Metabolism
"""


import odecell

import libsbml
import pandas as pd
import numpy as np

from collections import defaultdict

#########################################################################################
def initModel(sim_properties):
    """
    Initiate the model object and pass in the CME species counts

    Arguments:

    sim_properties['counts'] (particle map): The CME particle Map
    
    Returns:

    upModel (odecell model object): The updated ODE kinetic model for simulation
    """

    # Initialize the ODECell model
    model = odecell.modelbuilder.MetabolicModel()
    
    zeroOrderOnOff = '$onoff * $K'

    model.zeroOrderOnOff = odecell.modelbuilder.RateForm(zeroOrderOnOff)

    model.updateAvailableForms()

    # Set verbosity outputs to zero for now to improve performance
    model.setVerbosity(0)

    # Define Rxns and pass in the Particle Map containing enzyme concentrations
    model = defineRxns(model, sim_properties)

    return model
#########################################################################################


#########################################################################################
def defineRxns(model, sim_properties):
    
    model = addProteinMetabolites(model, sim_properties)
    
    model = defineRandomBindingRxns(model, sim_properties)
    
    model = defineNonRandomBindingRxns(model, sim_properties)
    
    # Other Random Bindig reactions are GLCK and GLCT reactions that converts extracellular glucose into g6p.
    # model = defineOtherRandomBindingReactions(model, sim_properties)
    
    return model
#########################################################################################


#########################################################################################
def partTomM(particles, sim_properties):
    """
    Convert particle counts to mM concentrations for the ODE Solver

    Parameters:
    particles (int): The number of particles for a given chemical species

    Returns:
    conc (float): The concentration of the chemical species in mM
    """

    ### Constants
    NA = 6.022e23 # Avogadro's

    conc = (particles*1000.0)/(NA*sim_properties['volume_L'][-1])

    return conc
#########################################################################################


#########################################################################################
def reptModel(model):
    """
    Report on the constructed hybrid model - but probably would only want to do after the first time step

    Arguments: 
    model (model obj.): The ODE Kinetic Model

    Returns:

    None
    """

    dictTypes = defaultdict(int)
    typeList = ["Transcription","Translation","Degradation"]

    for rxn in model.getRxnList():
        
        if rxn.getResult():
            # If an explicit result has been set, this is a dependent reaction.
            dictTypes["Dependent reactions"] += 1
            continue
        
        for rxntype in typeList:
            if rxntype in rxn.getID():
                dictTypes[rxntype] += 1

                
    #print( "There are {} ODEs in the model:".format(len(model.getRxnList())) )

    outList = list(dictTypes.items())
    outList.sort(key=lambda x: x[1], reverse=True)
    for key,val in outList:
        print("{:>20} :   {}".format(key,val) )
        return 0
    
    return None
#########################################################################################

# Central, Nucleotide, Lipid, Cofactor and Transprot Reactions were defined.
#########################################################################################
def defineRandomBindingRxns(model, sim_properties):
    """
    Define all of the reactions and rateforms needed for the current module to an existing module.

    """
    
    params_file = '../input_data/kinetic_params.xlsx'
    
    central_params = pd.read_excel(params_file, sheet_name='Central')
    nucleotide_params = pd.read_excel(params_file, sheet_name='Nucleotide')
    lipid_params = pd.read_excel(params_file, sheet_name='Lipid')
    cofactor_params = pd.read_excel(params_file, sheet_name='Cofactor')
    transport_params = pd.read_excel(params_file, sheet_name='Transport')
    
    metabolism_params = pd.concat([central_params, nucleotide_params, lipid_params, cofactor_params, transport_params], ignore_index=True) #, transport_params
    
    reaction_list = []

    for row, item in metabolism_params.iterrows():
        if item['Reaction Name'] not in reaction_list:
            reaction_list.append(item['Reaction Name'])


    # The .xml file for flux balance analysis
            
    sbmlFile = "../input_data/Syn3A_updated.xml"

    docSBML = libsbml.readSBMLFromFile(sbmlFile)
    modelSBML = docSBML.getModel()

    speciesNames = [spc.name for spc in modelSBML.getListOfSpecies()]
    speciesNamesLower = [x.lower() for x in speciesNames]
    speciesIDs = [spc.id for spc in modelSBML.getListOfSpecies()]

    rxnNamesSBML = [ x.name for x in modelSBML.getListOfReactions()]
    
    for rxnID in reaction_list:
    
#         print(rxnID)
        
        rxn_info = getSpecIDs(rxnID, modelSBML, rxnNamesSBML)
        
        # rxn_params are the rows with rxnID as reaction name
        rxn_params = metabolism_params.loc[ metabolism_params["Reaction Name"] == rxnID ]

#         print(rxn_info)

    #     print(rxn_params)

        # From .xml file to read the list of substrates and products and their stoichiometry
        substrates_list = rxn_info[0][0]
        substrates_stoich = rxn_info[0][1]
        products_list = rxn_info[1][0]
        products_stoich = rxn_info[1][1]
#         print(substrates_list)
#         print(substrates_stoich)
#         print(products_list)
#         print(products_stoich)
        
        substrate_count = int(-np.sum(substrates_stoich))
        product_count = int(np.sum(products_stoich))

        # Only need the number of reactions and products to determine the ratelaw
        rateLaw = Enzymatic(substrate_count, product_count)
        
#         print(rateLaw)
        
        rateName = rxnID+'_rate'
        
        # Do we need to addrateform per reaction?
        model.addRateForm(rateName, odecell.modelbuilder.RateForm(rateLaw))

        rxnIndx = model.addReaction(rxnID, rateName, rxnName="Reaction " + rxnID)

        kcatF = rxn_params.loc[ rxn_params["Parameter Type"] == "Substrate Catalytic Rate Constant" ]["Value"].values[0]
        kcatR = rxn_params.loc[ rxn_params["Parameter Type"] == "Product Catalytic Rate Constant" ]["Value"].values[0]
        
        model.addParameter(rxnIndx, 'kcatF', kcatF)
        model.addParameter(rxnIndx, 'kcatR', kcatR)
        
#         print(kcatF)
#         print(kcatR)

        rxn_KMs = rxn_params.loc[ rxn_params["Parameter Type"] == "Michaelis Menten Constant" ]
    
#         if rxnID == 'NADHK':
#             print(rxnID)
#             print(substrates_list)
#             print(substrates_stoich)
#             print(products_list)
#             print(products_stoich)
#             print(rateLaw)
        
        # Define the reactants, products, their concentrations, and Michaelis Menten Constants
        sub_rxn_indx_counter = 0
        for i in range(len(substrates_list)):

            metID = substrates_list[i]
            
#             if spcID.endswith("_e"):
            
            if metID not in list(model.getMetDict().keys()):
            
                if metID.endswith('_e'):
                    
                    # for excellular species, their concentration is fixed
                    spcConc = sim_properties['medium'][metID]
                    
                else:
                
                    spcConc = partTomM(sim_properties['counts'][metID][-1], sim_properties)
            
                    model.addMetabolite(metID, metID, spcConc)

            stoichiometry = int(-substrates_stoich[i])

            
            met_KM = rxn_KMs.loc[ rxn_KMs["Related Species"] == metID ]["Value"].values[0]
            
            # The stoichiometry of reactants can be more than 1
            for j in range(stoichiometry):
                
                sub_rxn_indx_counter = sub_rxn_indx_counter + 1
                
                rateFormID = 'Sub' + str(sub_rxn_indx_counter)
            
                if metID.endswith('_e'):
                    
#                     spcConc = sim_properties['medium'][metID]

#                     model.addParameter(rxnIndx, rateFormID, spcConc)
                    # for excellular species, their concentration is fixed and as parameter
                    model.addParameter(rxnIndx, rateFormID, sim_properties['medium'][metID])
#                     print(rxnIndx, rateFormID, metID)

                else:
#                     print(rxnIndx, rateFormID, metID)
                    if metID == 'M_o2_c':
                        model.addParameter(rxnIndx, rateFormID, metID)
                    else:
                        model.addSubstrate(rxnIndx, rateFormID, metID)
                
                KM_ID = 'KmSub' + str(sub_rxn_indx_counter)
                
                model.addParameter(rxnIndx, KM_ID, met_KM)
                
#                 print(rxnIndx, KM_ID, met_KM)




#             print(metID, stoichiometry, met_KM)
        # For product side
        prod_rxn_indx_counter = 0
        for i in range(len(products_list)):

            metID = products_list[i]
            
            if metID not in list(model.getMetDict().keys()):
                
                if metID.endswith('_e'):
                    
                    spcConc = sim_properties['medium'][metID]
                    
                else:
                
                    spcConc = partTomM(sim_properties['counts'][metID][-1], sim_properties)
                
                    model.addMetabolite(metID, metID, spcConc)
                
#             print('Added metabolite: ', metID, spcConc)

            stoichiometry = int(products_stoich[i])
#             print(metID, stoichiometry)

            met_KM = rxn_KMs.loc[ rxn_KMs["Related Species"] == metID ]["Value"].values[0]
            
            for j in range(stoichiometry):
                
                prod_rxn_indx_counter = prod_rxn_indx_counter + 1
                
                rateFormID = 'Prod' + str(prod_rxn_indx_counter)
            
                if metID.endswith('_e'):
                    
#                     spcConc = sim_properties['medium'][metID]

#                     model.addParameter(rxnIndx, rateFormID, spcConc)
                    model.addParameter(rxnIndx, rateFormID, sim_properties['medium'][metID])
#                     print(rxnIndx, rateFormID, metID)

                else:
#                     print(rxnIndx, rateFormID, metID)
                    model.addProduct(rxnIndx, rateFormID, metID)
#                     if j==0:
#                         if stoichiometry == 1:
#                             model.addProduct(rxnIndx, rateFormID, metID)
#                         else:
#                             print(rxnID, metID, stoichiometry)
#                             model.addProduct(rxnIndx, rateFormID, metID, stoich=int(stoichiometry))
#                     else:
#                         print(rxnID, metID, j)
#                         model.addParameter(rxnIndx, rateFormID, metID)
#                     print(j, rateFormID, metID)
                
                KM_ID = 'KmProd' + str(prod_rxn_indx_counter)
                
                model.addParameter(rxnIndx, KM_ID, met_KM)

#             print(metID, stoichiometry, met_KM)
            
        EnzymeConc = getEnzymeConc(rxn_params, sim_properties)
            
#         EnzymeConc = partTomM(rxn_params.loc[ rxn_params["Parameter Type"] == "Eff Enzyme Count" ]["Value"].values[0], sim_properties)
        
        model.addParameter(rxnIndx, "Enzyme", EnzymeConc)
        
        model.addParameter(rxnIndx, "onoff", 1, lb=0, ub=1)
        
#         print(' ')
            

#     print("Reactions defined")

    reptModel(model)

    return model
#########################################################################################


#########################################################################################
def defineOtherRandomBindingReactions(model, sim_properties):
    # For other random binding reactions, we give the information of metID and stoichiometry in the excel file already
    params_file = '../input_data/kinetic_params.xlsx'
    
    RXNS_params = pd.read_excel(params_file, sheet_name='Other-Random-Binding')
    
    reaction_list = []

    for row, item in RXNS_params.iterrows():
        if item['Reaction Name'] not in reaction_list:
            reaction_list.append(item['Reaction Name'])
#             print(item['Reaction Name'])
            
    for rxnID in reaction_list:
    
#         print(rxnID)

        rxn_params = RXNS_params.loc[ RXNS_params["Reaction Name"] == rxnID ]
        
        substrate_count = int(rxn_params.loc[ rxn_params["Parameter Type"] == "Substrates" ]["Value"].values[0])
        product_count = int(rxn_params.loc[ rxn_params["Parameter Type"] == "Products" ]["Value"].values[0])
        
        rateLaw = Enzymatic(substrate_count, product_count)
        
        rateName = rxnID+'_rate'
        
        model.addRateForm(rateName, odecell.modelbuilder.RateForm(rateLaw))

        rxnIndx = model.addReaction(rxnID, rateName, rxnName="Reaction " + rxnID)

        kcatF = rxn_params.loc[ rxn_params["Parameter Type"] == "Substrate Catalytic Rate Constant" ]["Value"].values[0]
        kcatR = rxn_params.loc[ rxn_params["Parameter Type"] == "Product Catalytic Rate Constant" ]["Value"].values[0]
        
        model.addParameter(rxnIndx, 'kcatF', kcatF)
        model.addParameter(rxnIndx, 'kcatR', kcatR)
        
#         print(kcatF)
#         print(kcatR)

        rxn_KMs = rxn_params.loc[ rxn_params["Parameter Type"] == "Michaelis Menten Constant" ]
        
        
        for i in range(1, substrate_count+1):

            metID = rxn_params.loc[ rxn_params["Parameter Type"] == "Sub" + str(i) ]["Value"].values[0]
            
#             if spcID.endswith("_e"):
            
            if metID not in list(model.getMetDict().keys()):
            
                if metID.endswith('_e'):
                    
                    spcConc = sim_properties['medium'][metID]
                    
                else:
                
                    spcConc = partTomM(sim_properties['counts'][metID][-1], sim_properties)
                
#                 print('Added metabolite: ', metID, spcConc)
                
                    model.addMetabolite(metID, metID, spcConc)

            stoichiometry = int(1)
            
            met_KM = rxn_KMs.loc[ rxn_KMs["Related Species"] == metID ]["Value"].values[0]
            
#             for j in range(stoichiometry):
                
#                 sub_rxn_indx_counter = sub_rxn_indx_counter + 1
                
            rateFormID = 'Sub' + str(i)
            
            if metID.endswith('_e'):

#                     spcConc = sim_properties['medium'][metID]

#                     model.addParameter(rxnIndx, rateFormID, spcConc)
                model.addParameter(rxnIndx, rateFormID, sim_properties['medium'][metID])
#                     print(rxnIndx, rateFormID, metID)

            else:
#                     print(rxnIndx, rateFormID, metID)
                model.addSubstrate(rxnIndx, rateFormID, metID)

            KM_ID = 'KmSub' + str(i)

            model.addParameter(rxnIndx, KM_ID, met_KM)
                

        for i in range(1, product_count+1):

            metID = rxn_params.loc[ rxn_params["Parameter Type"] == "Prod" + str(i) ]["Value"].values[0]
            
#             if spcID.endswith("_e"):
            
            if metID not in list(model.getMetDict().keys()):
            
                if metID.endswith('_e'):
                    
                    spcConc = sim_properties['medium'][metID]
                    
                else:
                
                    spcConc = partTomM(sim_properties['counts'][metID][-1], sim_properties)
                
#                 print('Added metabolite: ', metID, spcConc)
                
                    model.addMetabolite(metID, metID, spcConc)

            stoichiometry = int(1)
            
            met_KM = rxn_KMs.loc[ rxn_KMs["Related Species"] == metID ]["Value"].values[0]
            
#             for j in range(stoichiometry):
                
#                 sub_rxn_indx_counter = sub_rxn_indx_counter + 1
                
            rateFormID = 'Prod' + str(i)
            
            if metID.endswith('_e'):

#                     spcConc = sim_properties['medium'][metID]

#                     model.addParameter(rxnIndx, rateFormID, spcConc)
                model.addParameter(rxnIndx, rateFormID, sim_properties['medium'][metID])
#                     print(rxnIndx, rateFormID, metID)

            else:
#                     print(rxnIndx, rateFormID, metID)
                model.addProduct(rxnIndx, rateFormID, metID)

            KM_ID = 'KmProd' + str(i)

            model.addParameter(rxnIndx, KM_ID, met_KM)

            
        EnzymeConc = getEnzymeConc(rxn_params, sim_properties)
            
#         EnzymeConc = partTomM(rxn_params.loc[ rxn_params["Parameter Type"] == "Eff Enzyme Count" ]["Value"].values[0], sim_properties)
        
        model.addParameter(rxnIndx, "Enzyme", EnzymeConc)
        
        model.addParameter(rxnIndx, "onoff", 1, lb=0, ub=1)
        
    return model

#########################################################################################


#########################################################################################
def defineNonRandomBindingRxns(model, sim_properties):
    # For non random binding reactions, everything is given in the excel sheet
    params_file = '../input_data/kinetic_params.xlsx'
    
    RXNS_params = pd.read_excel(params_file, sheet_name='Non-Random-Binding Reactions')
    
    reaction_list = []
    # Get the full list of reaction names
    for row, item in RXNS_params.iterrows():
        if item['Reaction Name'] not in reaction_list:
            reaction_list.append(item['Reaction Name'])
#             print(item['Reaction Name'])
            
    for rxnID in reaction_list:
    
#         print(rxnID)

        rxn_params = RXNS_params.loc[ RXNS_params["Reaction Name"] == rxnID ]
        # Ratelaws are given mannually
        rateLaw = str(rxn_params.loc[ rxn_params["Parameter Type"] == "Kinetic Law" ]["Value"].values[0])
        
        rateName = rxnID+'_rate'
        
        model.addRateForm(rateName, odecell.modelbuilder.RateForm(rateLaw))

        rxnIndx = model.addReaction(rxnID, rateName, rxnName="Reaction " + rxnID)
        
        for index, row in rxn_params.iterrows():
            
            param = row['Parameter Type']
            
            if (param != "Reaction Formula") and (param != "Kinetic Law"):

                if param.startswith('Sub'):
                    
                    metID = row['Value']
                    
                    if metID not in list(model.getMetDict().keys()):
            
                        if metID.endswith('_e'):

                            spcConc = sim_properties['medium'][metID]

                        else:

                            spcConc = partTomM(sim_properties['counts'][metID][-1], sim_properties)

#                         print('Added metabolite: ', metID, spcConc)

                            model.addMetabolite(metID, metID, spcConc)

                    if metID.endswith('_e'):

                        model.addParameter(rxnIndx, param, sim_properties['medium'][metID])

                    else:

                        model.addSubstrate(rxnIndx, param, metID)

                elif param.startswith('Prod'):
                    
                    metID = row['Value']
                    
                    if metID not in list(model.getMetDict().keys()):
            
                        if metID.endswith('_e'):

                            spcConc = sim_properties['medium'][metID]

                        else:

                            spcConc = partTomM(sim_properties['counts'][metID][-1], sim_properties)

#                         print('Added metabolite: ', metID, spcConc)

                            model.addMetabolite(metID, metID, spcConc)

                    if metID.endswith('_e'):

                        model.addParameter(rxnIndx, param, sim_properties['medium'][metID])

                    else:

                        model.addProduct(rxnIndx, param, metID)

                else:

                    model.addParameter(rxnIndx, param, row['Value'])

    return model
#########################################################################################

# Protein not as enzymes but as reactants and products; Phosphorelay: 4 proteins; P_0227 PDH
# BaseForm in the original form of protein, listed as the first one in Metabolite IDs
# This function only defines the names of metabolites; the reactions and rates are included in others.
#########################################################################################
def addProteinMetabolites(model, sim_properties):
    """
    
    Description: add the counts of each form of proteins to the ODE model
    
    For P_0621, P_0065, P_0227, their initial counts are given by initializeMetabolitesCounts
    For 4 proteins in phosphorelay, their initial counts are given by initializeProteinMetabolitesCounts

    The new generated protein from CME side will add to the base form ptsi ... in ODE
    """ 


    data_file = '../input_data/initial_concentrations.xlsx'
    
    ptnMets = pd.read_excel(data_file, sheet_name='protein_metabolites')
    
    for index, row in ptnMets.iterrows():
        
        ptnID = row['Protein']
        
        metabolites = row['Metabolite IDs'].split(',')
        

        # updating the total protein counts from CME side per second
        ptnCount = sim_properties['counts'][ptnID][-1]
        
        formsCount = 0
        
        for metID in metabolites:
            
            if metID != metabolites[0]:
                
                formsCount = formsCount + sim_properties['counts'][metID][-1]
                
#                 print(metID, formsCount)
                
                model.addMetabolite(metID, metID, partTomM(sim_properties['counts'][metID][-1], sim_properties))
                
        baseFormCount = int(ptnCount - formsCount)
        
#         print(metabolites[0], baseFormCount)
        
        baseFormConc = partTomM(baseFormCount, sim_properties)

        model.addMetabolite(metabolites[0], metabolites[0], baseFormConc)
        
    return model
#########################################################################################

# Enzymatic: Define the form of Random Binding Rate Constant
#########################################################################################
def Enzymatic(subs, prods):
        
    def numerator(subs,prods):
        
        subterm = [ '( $Sub' + str(i) + ' / $KmSub' + str(i) + ' )' for i in range(1,subs+1)]
        subNumer = ' * '.join(subterm)
        
        prodterm = [ '( $Prod' + str(i) + ' / $KmProd' + str(i) + ' )' for i in range(1,prods+1)]
        prodNumer = ' * '.join(prodterm)
        
        numerator = '( ' + '$kcatF * ' + subNumer + ' - ' + '$kcatR * ' + prodNumer + ' )'
        return numerator
    
    def denominator(subs,prods):
        
        subterm = [ '( 1 + $Sub' + str(i) + ' / $KmSub' + str(i) + ' )' for i in range(1,subs+1)]
        subDenom = ' * '.join(subterm)
        
        prodterm = [ '( 1 + $Prod' + str(i) + ' / $KmProd' + str(i) + ' )' for i in range(1,prods+1)]
        prodDenom = ' * '.join(prodterm)
        
        denominator = '( ' + subDenom + ' + ' + prodDenom + ' - 1 )'
        return denominator
        
    rate = '$onoff * $Enzyme * ( ' + numerator(subs,prods) + ' / ' + denominator(subs,prods) + ' )'
    
    return rate
#########################################################################################


#########################################################################################
def getEnzymeConc(rxn_params, sim_properties):
    
    EnzymeStr = rxn_params.loc[ rxn_params["Parameter Type"] == "Eff Enzyme Count" ]["Value"].values[0]
    
    Enzymes = EnzymeStr.split('-')
#     print(Enzymes)
    
    if len(Enzymes) == 1:
        
        if Enzymes[0] == 'default':
            
#             EnzymeConc = partTomM(10, sim_properties)
            EnzymeConc = 0.001
            
            return EnzymeConc
        
        else:
            
            ptnID = Enzymes[0]
            
#             print(sim_properties['counts'][ptnID])
            
            EnzymeConc = partTomM(sim_properties['counts'][ptnID][-1], sim_properties)
            
            return EnzymeConc
    
    else:
        # Here we difine two means of how multiple enzymes works
        # For 'or' case, all enzymes can catalyze the reaction and we give them the same kinetic parameter
        # And the rate is determined by the sum of the concentrations of enzymes
        GPRrule = rxn_params.loc[ rxn_params["Parameter Type"] == "GPR rule" ]["Value"].values[0]

        if GPRrule == 'or':
            
            ptnCount = 0
            
            for ptnID in Enzymes:
                
                ptnCount = ptnCount + sim_properties['counts'][ptnID][-1]
                
            EnzymeConc = partTomM(ptnCount, sim_properties)
                
            return EnzymeConc

        # For 'and' case, all enzymes are needed to catalyze the reaction and the concentration of the lowest enzyme determines the rate        
        elif GPRrule == 'and':
            
            ptnCounts = []
            
            for ptnID in Enzymes:
                
                ptnCounts.append(sim_properties['counts'][ptnID][-1])
                
            ptnCount = min(ptnCounts)
            
            EnzymeConc = partTomM(ptnCount, sim_properties)
                
            return EnzymeConc
        
    print('Something went wrong getting enzyme count')
        
    return None
#########################################################################################

# From .xml file to obtain the metID and stoichiometry information
#########################################################################################
def getSpecIDs(rxnName, modelSBML, rxnNamesSBML):
    
    returnList = []
    
    rxnObj = modelSBML.getReaction( rxnNamesSBML.index(rxnName) )
    
    # Use model SBML to get IDs, names, and stoichiometries for reactants
    specIDs = [ x.getSpecies() for x in rxnObj.getListOfReactants() ]
    spcStoich = [ -1*float(x.getStoichiometry()) for x in rxnObj.getListOfReactants() ]
    spcNames = [ modelSBML.getSpecies( spcID ).name for spcID in specIDs]
    
    specIDs_noH = []
    spcStoich_noH = []
    
    for i in range(len(specIDs)):
        
        metID = specIDs[i]
        
        if (metID != 'M_h_c') and (metID != 'M_h_e') and (metID != 'M_h2o_c')  and (metID != 'M_h2o_e'):
            
            specIDs_noH.append(metID)
            stoich = spcStoich[i]
            spcStoich_noH.append(stoich)
    
    if np.any( np.isnan( spcStoich ) ):
        raise Exception('Invalid stoichiometry for reaction: {}'.format(rxnName)) 
    
#     returnList.append( [specIDs, spcStoich] )
#     returnList.append( [spcNames, specIDs, spcStoich] )
    returnList.append( [specIDs_noH, spcStoich_noH] )
    
    # Now do the same for products
    specIDs = [ x.getSpecies() for x in rxnObj.getListOfProducts() ]
    spcStoich = [ float(x.getStoichiometry()) for x in rxnObj.getListOfProducts() ]
    spcNames = [ modelSBML.getSpecies( spcID ).name for spcID in specIDs]
    
    specIDs_noH = []
    spcStoich_noH = []
    
    for i in range(len(specIDs)):
        
        metID = specIDs[i]
        
        if (metID != 'M_h_c') and (metID != 'M_h_e') and (metID != 'M_h2o_c')  and (metID != 'M_h2o_e'):
            
            specIDs_noH.append(metID)
            stoich = spcStoich[i]
            spcStoich_noH.append(stoich)
    
    if np.any( np.isnan( spcStoich ) ):
        raise Exception('Invalid stoichiometry for reaction: {}'.format(rxnName)) 
    
#     returnList.append( [specIDs, spcStoich] )
#     returnList.append( [spcNames, specIDs, spcStoich] )
    returnList.append( [specIDs_noH, spcStoich_noH] )
    
    return returnList
#########################################################################################


