### syn3A.gb

[Genbank file](https://www.ncbi.nlm.nih.gov/nuccore/CP016816) of JCVI-SYN3A (ACCESSION: CP016816 in NCBI)

syn3A.gb contains the sequence, segmentation, functions of the whole genome of syn3A.

Read in once at the very beginning of the simulation

**mapDNA** function in **initiation.py** will convert syn3A.gb into a multilayer dictionary *genome*, which will be iterated to define the CME trascription, translation, and degradation reactions in **addGeneticInformationReactions** function in **rxns_CME.py**.

### Syn3A_updated.xml

Updated SBML file of syn3A based on the file from [elife 2019 paper](https://elifesciences.org/articles/36842).

Syn3A_updated.xml contains compartments, metabolties, reactions, and related genes of metabolism in JCVI-Syn3A.

Read in once when constructing ODEs in **rxns_ODE.py** for ***only the stoichiometric coefficients*** of reactions.

### initial_concentration.xlsx

Excel file contains the counts/concentrations of proteins, medium and metabolites.

Read in once when initialize the counts/concentrations of species at the beginning of the simulation 

Different sheets:

1) **Comparative Proteomics**: 

    Read in when adding initial counts of proteins in CME

    Gene LocusTag and the function, localization, and initial counts of its protein

2) **Experimental Medium**:

    Not Read in

    The concentrations experimental medium

3) **Simulation Medium**:

    Read in when initializing the medium

    The concentrations, Met ID of simulated medium

4) **Intracellular Metabolites**

    Read in when initializing the cytoplasmic metabolties

    The concentration, Met ID of simulated metabolites

5) **protein_metabolites**

    Read in when initializing protein metabolties in **initiation.py** and when adding counts of protein metabolits in **rxns_ODE.py**

    The Met IDs of proteins

### kinetic_params.xlsx

Excel file contains all the reactions in ODE and tRNA charging in CME.

Read in when constructing the reactions in ODEs (per hookInterval) and tRNA charging in CME

**Central**, **Nucleotide**, **Lipid**, **Cofactor**, and **Transport** sheets contain reactions obeying random binding model and convenience rate law. These reactions will be added via **defineRandomBindingRxns** function in **rxns_ODE.py**.

**Non-Random-binding Reactions** sheet contains literally non-random-binding reactions. They are serial phosphorelay reactions and passive transport reactins. These reactions will be added via **defineNonRandomBindingRxns** function in **rxns_ODE.py**.

**tRNA Charging** sheet contains amino acid, synthetase, and kinetic parameters for tRNA charging reactions and will be read in via **tRNAcharging** function in rxns_CME.py

Other sheets in **kinetic_params.xlsx** are not used in current simulation.

