## Programs for running CME/ODE model of minimal cell

#### Author: Enguang Fu
#### Email: enguang3@illinois.edu


#### Generals      

Simulating a whole Syn3A cell using hybrid CME/ODE method

The CME and ODE communicates per hookInterval using user defined solver hookSolver_CMEODE

Launch the simulation by mpirun.sh



#### Input

Folder: ../input_data/

Check Input_README.md for more

syn3A.gb

    gene bank file of Syn3A; 
    converted into multilayer dictionary containing all genetic information by mapDNA function in initiation.py

initial_concentrations.xlsx

    initials counts of proteins; concentrations of metabolits in CME
    called when initializing species in CME and ODE

kinetic_params.xlsx

    kinetic parameters for metabolic reactions in ODE and tRNA charging in CME
    called by rxns_CME.py and rxns_ODE.py to get the reactions' parameters 


Syn3A_updated_xml

    Systems Biology Markup Language (SBML) file of Syn3A 
    called by rxns_ODE.py to get stoichiometric coefficient of metabolic reactions in ODE

#### Scripts

Main Script:

    WCM_CMEODE_Hook.py - Main script

Launching the Parallel Simulation:

    mpirun.sh - simulation launching bash file

Constructing the CMEODE framework:

    species_counts.py - defining class for easy passing of species counts data in the hook algorithm.

    integrate.py - perform the ODE integration 

    initiation.py - initialize constants, time traces of counts, membranes 

    communicate.py - updates CME and ODE states, calculate costs and update membrane
    
    hookSolver_CMEODE.py - user defined solver where hookSimulation performs every hookInterval

    hook_CMEODE.py - the operations when do hooking/communicating between CME/ODE

    filesaving.py - export the time traces of counts, Surface area, and fluxes into CSV files

Biological:

    GIP_rates.py - calculate the rates for genetic information processes (GIP) in CME

    rxns_ODE.py - construct the ODE system using odecell

    rxns_CME.py - add GIP reactions including replication, transcription, translation, tRNA charging to CME

    replication_multiple.py - replication initiation and replication reactions



#### Output

Output folders are defined in the mpirun.sh and will be created if not exist.

The output files are three csv files and one log file for each CMEODE simulation replicate with rank r.

    Three csv files are counts_r.csv, SA_r.csv, and Flux_r.csv and one log file is log_r.txt.

    counts_r.csv is the trajectories of counts (unit: numbers) of all CMEODE species.

    SA_r.csv is the trajectories of surface area (unit: nm2 or m2) and volume (unit: Liter).

    Flux_r.csv is the trajectorie of all fluxes (unit: mM/s) through ODE reactions.
 
    The trajectories are sampled per hook interval.

    The typical size for a single csv file is 10 to 200 MBs when simulating 6300 seconds with hook interval 1 second.

log_r.txt as a log file contains the simulation start and end time, time cost of each operations and possible errors.




#### Launching the Simulation


Execute the mpirun.sh to parallelly run CMEODE models 

See mpi_run.sh file for detailed explanation


#### Python Environment

Lattice Microbe

odecell

    construct the ODEs of the metabolism

mpi4py

    run parallel simulations



#### Explanation


genome = sim_properties['genome']
genome: 
Key: JCVISYN3A_****; Total number is 496
Value: A subdictionary with keys Type, startIndex, endIndexm originalEnd, RNAsquence, AAsequence, GeneName in JCVISYN3A_****

One entry in genome dictionary
JCVISYN3A_0250': {'Type': 'protein',
  'startIndex': [42472],
  'originalStart': 42472,
  'endIndex': [42403],
  'originalEnd': 42403,
  'RNAsequence': 'AUGAAUAAAAAAGAGAUAUUUAACACUGAUUUUUUUGAAUCAGGUUUAGCUUAUAUUCUAACUAAUUUAGAUUUUAUUCAAGAAGAAUUAGAACAAGAAAAACUACAAACUAGUCUAGUAGAAAAAUU
AAUAACUGAUUUUGAAGAUGUUGAAGAUUAUGAAACAUGAGAUUUAUUAACUAAUAAUUUAAUUCAAUCUGAAGAUAAGAUUUUAGAAGAAAUUCAAAAAAUAAAAGACUCAACAAAAUUCAAUUUAUUAAAUAGUUAUUUUUUAGCA
AAAAACCUUGCUAUUUAUUUAAAGUCUAAUAGUUUUUUAAUAGAACAAAUAAACAAAUUACAAACAAAUUCUCCUGAUGAUUUAUCAGAAGAUAAAAAAGAAGAAUUUAUUAAUAACUUAAAACAAGAAAUUUUAAAAAAUAAUUCAGA
AUUAUACAAACAAAACGAAAGAUUAUUUAAAGAAAUAUUUGAUAAAAAAGUUGAAUUUAAAAAAAUCUAUCAACUUUUAAUUAAAGAAACUGAGUUUGAAGAUUUUAAUUAUGCAAACGAAUUAUUAUUUAAUAUGUUAAAUAACAACUUCA
AAUUUAAUAACAAACAAGAUUUAUUAAAAUUAGAAGUUUUAAAUAAUGCUCAAUCUUUAAUAGAUUUUCUUACUUUUUAUGAAUCUAGUUUAUUUGAUGAUGAAAAAGAAUAA',
  'AAsequence': 'MNKKEIFNTDFFESGLAYILTNLDFIQEELEQEKLQTSLVEKLITDFEDVEDYETWDLLTNNLIQSEDKILEEIQKIKDSTKFNLLNSYFLAKNLAIYLKSNSFLIEQINKLQTNSPDDLSEDKKEEFINNLKQE
ILKNNSELYKQNERLFKEIFDKKVEFKKIYQLLIKETEFEDFNYANELLFNMLNNNFKFNNKQDLLKLEVLNNAQSLIDFLTFYESSLFDDEKE*',
  'GeneName': 'Uncharacterized protein'}

Type = ['protein', 'ncRNA', 'gene', 'rRNA', 'tRNA', 'tmRNA']
455 proteins, 2 ncRNA, 3 pseudo genes, 6 rRNA, 29 tRNA, 1 tmRNA


Recording the Log file:
Many solution. From bash side, from MPI side, from the working python script side
Here we adopt the way that redirecting the sys.stdout and sys.stderr in one pytonto one log file 
