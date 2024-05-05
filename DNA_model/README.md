# Hands-on tutorial: Simulating DNA replication and dynamics with btree_chromo and LAMMPS
## Description
This tutorial was prepared for Day 1 of the STC QCB Advanced Computational Workshop, held on Monday May 6 2024.

Here, we simulate DNA replication and dynamics using the program btree_chromo, available online at  [github.com/brg4/btree_chromo](https://github.com/brg4/btree_chromo).

Powerpoint slides to go with tutorial instructions:
[github.com/Luthey-Schulten-Lab/Workshop_2024/tree/main/DNA_model](https://github.com/Luthey-Schulten-Lab/Workshop_2024/tree/main/DNA_model)

**Outline of tutorial:** \
Part 1: Modeling Replication States \
Part 2: Preparing the Physical Structure \
Part 3: Simulating Chromosome Dynamics \
Part 4: Analysis, Conversion to LM and MARTINI (slides only)

## Setting up the tutorial
**Step 1: Log in to Delta**
```bash
ssh $USERNAME@login.delta.ncsa.illinois.edu
```
Replace the $USERNAME with your username. You will need to type your password and do 2FA.

**Step 2: Create workspace and copy examples folder**
```bash
bash /projects/bcuj/sharefile/Workshop_2024/DNA_model/prelaunch_btree_chromo.sh  
```
This bash script creates a workspace `/projects/bcuj/${USER}/btree_chromo_workspace`. It also copies the `examples` directory into the workspace, which contains example input and output files for **btree_chromo**. 

**Step 3: Launch the container**
```bash
bash /projects/bcuj/sharefile/Workshop_2024/DNA_model/launch_btree_chromo.sh  
```
This will run the container in interactive mode. You should now see the `Apptainer>` prompt which indicates your have entered the container. We will be running btree_chromo and viewing the terminal output in the container.

**Important:** The btree_chromo  exectuable is within the container in `/Software/btree_chromo/build/apps/`. \
**Important:** We have mounted `/projects/bcuj/${USER}/btree_chromo_workspace/examples` into the container in `/mnt/examples`. All changes made in the workspace will be reflected in the container.

**Step 4: Open a new terminal window**\
In a new terminal window, repeat step 1 and `cd  /projects/bcuj/${USER}/btree_chromo_workspace/examples`. Viewing and editing of the example files will be done within this terminal window.
## Modeling Replication States
In this section, you will learn how to represent replication states, including nested theta structures, with a binary tree model.

In your second terminal, `ls /preparing_chromosome`:
|File name| Description |
| --- | --- |
| preparing_chromosome_directives.inp|contains directives (commands to be executed by btree_chromo) |
| chromo_state_0.dat | contains the initial replication state |
| chromo_state_1.dat | contains the final replication state |
| transforms.dat | contains a set of replication events |

Using a text exitor such as **vim**, open 'preparing_chromosome_directives.inp'. Scroll through the file to see all of the commands. Here, **btree_chromo** first creates a chromosome with 1000 monomers, and then applies a series of transforms, printing the replication state at each step, and finally outputs the final state.

Make sure that the "terminate'' at the top is commented to execute all of the commands.

Next, in the Apptainer window:
```bash
cd /Software/btree_chromo/build/apps
```
and run 'preparing_chromosome_directives.inp':
```bash
./btree_chromo /mnt/examples/preparing_chromosome/preparing_chromosome_directives.inp
```
Check out the terminal output. Notice how **btree_chromo** starts by testing the validity of all the commands and parameters and prints them to the terminal output. It then executes each of the commands. For example, you should see this output corresponding to the first replication transformation:

```bash
COMMAND: transform
	param_0: m_cw100_ccw200



COMMAND: print



printing tree with 2 leaves and 1 forks
fork breakdown: 0 completed, 1 active
leaves: 
ml mr 
active forks: 
m 
total_size = 1300
| generation = 0
| rho_t = 300/1000, rho_cw = 100/800, rho_ccw = 200/900
| start = 0, mid = 500, end = 999
| start_link = 999, end_link = 0
| left branch
  | generation = 1
  | start = 0, mid = 500, end = 999
  | start_link = 999, end_link = 0
| right branch
  | generation = 1
  | start = 1000, mid = 1200, end = 1299
  | start_link = 299, end_link = 600
```


**Questions** 
1. Comment out all the transforms except for the one on line 76, which reads in the file 'transforms.dat'. Exit the directives file and open 'transforms.dat'. Based on these transforms, **what should the final size of the chromosome be?**   Run ‘preparing_chromosome_directives.inp’ and check the terminal output for the 'print' command to see if you were right. 
2. Now, apply your understanding: edit 'transforms.dat' to **create a scenario where the DNA is exactly doubled**. Challenge: How many leaves and generations can you make?

## Preparing the Physical Structure
In this section, you will prepare input coordinates for monomers confined by boundary particles and avoiding ribosomes, and will obtain daughter monomer positions using the train track model. (No actual dynamics will be simulated in this section.)

In your second terminal, move to the next example in `examples/preparing_physical_structure`.
There are several file types in this directory. Here are what they are used for:
|File type| Description |
| --- | --- |
| .inp| directives (commands to be executed by btree_chromo) |
| .dat | replication state |
| .xyz | monomer coordinates for quick visualization |
| .bin | monomer/ribosome coordinates (binary) |
| data. | LAMMPS data file |

Edit line 62 in 'preparing_physical_structure_directives.inp' so that it uses your transforms file that you made in the previous section. It should now look like this:
```bash
transforms_file:/mnt/examples/preparing_chromosome/transforms.dat
```
Make sure you are saving your changes.

In the Apptainer window run 'preparing_physical_structure_directives.inp':
```bash
./btree_chromo /mnt/examples/preparing_chromosome/preparing_physical_structure_directives.inp
```
**Visualizing train-track replication with VMD:**\
You will now copy over the .xyz files from Delta to your local machine in order to visualize them in vmd.
```bash
scp $USERNAME@login.delta.ncsa.illinois.edu:/projects/bcuj/$USERNAME/btree_chromo_workspace/examples/preparing_physical_structure/*.xyz .
```
(Alternatively, instead of `.` you may choose to specify path to a local directory.)

Open vmd and load a new molecule with the unreplicated coordinates (File->New Molecule...). Then load a new molecule with the replicated coordinates. Use the D (Draw) toggle in the VMD Main Window to see the replicated monomers appear and dissapear. No dynamics have been performed, so they should be close to their mothers!

Feel free to change the representation of the DNA monomers. We recommend using VDW representation style and increasing their radii to something bigger like 13. However, since the main point of this section was a quick visualization, don't worry about it too much. In the next section, we will provide a .tcl script which will create nice representations for the DNA, ribosomes and boundary particles.  

## Simulating Chromosome Dynamics
Take a look in the provided README for btree_chromo: [github.com/brg4/btree_chromo](https://github.com/brg4/btree_chromo). Scroll past the description and installation steps to see the list of possible directives.
Let's focus our attention on the commands that we can use to toggle various aspects of the Brownian dynamics simulation in LAMMPS, under Spatial System for Simulations: bonds, bending angles, and twisting angles.
| Directive | Description |
| --- | --- |
| `switch_bonds:(T/F)` | enable/disable bonds between DNA monomers (default T)
 | `switch_bending_angles:(T/F)` | enable/disable bending angles between DNA monomers (default T)
 | `switch_twisting_angles:(T/F)` | enable/disable twisting angles between DNA monomers (default T)

*Note: these directives must be used prior to 'write_LAMMPS_data_file' to take effect*.
By turning off these switches, we are effectively removing terms from the whole energy function that correspond to adjacent-monomer interactions in the DNA polymer: 
$$U= \sum_{i=1}^{N_{\mathrm{DNA}}}\left[U_i^b+U_i^t+U_i^a+U_i^s\right] +\sum_{i=1}^{N_{\mathrm{DNA}}-1} \sum_{j=i+1}^{N_{\mathrm{DNA}}} U_{i j}^{\mathrm{DNA}-\mathrm{DNA}}+\sum_{i=1}^{N_{\mathrm{DNA}}} \sum_j^{N_{\text {ribo }}} U_{i j}^{\mathrm{DNA}-\text { ribo }} +\sum_{i=1}^{N_{\text {ribo }}-1} \sum_{j=i+1}^{N_{\text {ribo }}} U_{i j}^{\text {ribo-ribo }} +\sum_{i=1}^{N_{\text {bdry }}} \sum_j^{N_{\mathrm{DNA}}} U_{i j}^{\text {bdry-DNA }}+\sum_{i=1}^{N_{\text {bdry }}} \sum_j^{N_{\text {ribo }}} U_{i j}^{\text {bdry-ribo }}.$$
Turning off bending removes the $U_i^b$ (cosine potential for bending), turning off twisting removes the $U_i^t$ and $U_i^a$ (cosine potentials for twisting and aligning), and turning off bonds will remove all of these as well as the $U_i^s$ (FENE potentials for stretching) resulting in separate diffusing DNA monomers rather than a DNA polymer. The other five terms in the energy function are for excluded volume interactions (purely repulsive Weeks-Chandler-Andersen (WCA) pair potentials).

We will now perform several simulation runs and will visualize their trajectories in VMD: 
1. Simulation with replication, using full Hamiltonian (the energy function with all terms) 
2. Simulation without twisting and bending potentials 
3. Simulation without bonds to show diffusing DNA monomers
   
Each run should not take more than 3 minutes.

**Full Hamiltonian run:**\
In your second terminal, go to `/examples/simulating_chromosome_with_replication/` and open 'simulating_chromosome_with_directives.inp'. Make sure that "terminate" is commented, and "switch_skip_runs" is set to "F". 

In line 48, rename the LAMMPS simulation output to ",example_full".

*Optional: Edit line 87 in 'simulating_chromosome_with_directives.inp' so that it uses your transforms file that you made in the previous section. It should now look like this:* `transforms_file:/mnt/examples/preparing_chromosome/transforms.dat`. *If you decide to do this, you may have to change line 78 to* `repeat:1`.

Save your changes and run the directives file. 

**No twist/no bend run:**\
In line 48, rename the LAMMPS simulation output to ",example_flexible".
Below line 48, add the lines "switch_bending_angles:F", and "switch_twisting_angles:F".

Save your changes and run the directives file. 

**Unbonded run:**\
In line 48, rename the LAMMPS simulation output to ",example_unbonded".
Below line 48, add the line "switch_bonds:F".

Save your changes and run the directives file. 

**Visualizing LAMMPS trajectories with VMD:**\
You will now copy over the .lammpstrj files from Delta to your local machine in order to visualize them in vmd.
```bash
scp $USERNAME@login.delta.ncsa.illinois.edu:/projects/bcuj/$USERNAME/btree_chromo_workspace/examples/simulating_chromosome_with_replication/*.lammpstrj .
```
(Alternatively, instead of `.` you may choose to specify path to a local directory.)
We will also copy over the .tcl scripts which will create nice representations for the DNA, ribosomes and boundary particles.  
```bash
scp $USERNAME@login.delta.ncsa.illinois.edu:/projects/bcuj/sharefile/Workshop_2024/DNA_model/*.tcl .
```

In VMD, delete the previous two molecules. Open the VMD TkConsole and do (Extensions->Tk Console). In the Tk Console, do `source load_example_full.tcl`. See the DNA polymer replicate!
Then do `source load_example_flexible.tcl` and `source load_example_unbonded.tcl` and note the differences.
|Monomer type|Color|Size|
|---|---|---|
|DNA|viridis color scale|12.0|
|_Ori_|red|20.0|
|_Ter_|orange|20.0|
|Fork|magenta|20.0|
|Ribosome|mauve|20.0|
|Boundary|gray|20.0|

## Simulating Chromosome Dynamics: SMC looping and topoisomerase

Finally, let's run simulations of an unreplicated chromosome under the influence of loops and topoisomerase. The relevant examples are in `/examples/simulating_with_loops_and_topo`. The relevant directives are described in the [README](https://github.com/brg4/btree_chromo/) under Simulator:
| Directive | Description |
| --- | --- |
| `simulator_load_loop_params:loop_params_file` | read a file (loop_params_file) containing the parameters for the looping interactions|
|`simulator_run_loops:Nloops,Nsteps,Tfreq,Dfreq,append_option,skip_option`|run Brownian dynamics with the hard/FENE potential for Nsteps with Nloops randomly placed, while printing thermodynamic information every Tfreq steps and dumping every Dfreq steps - append_option = noappend/append and skip_option = first/skip_first|

If we take a look at loop_params.txt, we find that it specifies the minimum number of monomers separating anchor and hinge, the distribution of extrusion steps, hinge unbinding probability and grab radius, loop update frequency and topoisomerase relaxation frequency.

```bash
# system parameters

# minimum distance separating anchor and hinge on strand - [# monomers]
min_dist=5

# distribution family for steps
family=poisson
# average extrusion distance [# monomers]
ext_avg=20.0
# max extrusion distance [# monomers]
ext_max=30

# hinge unbinding probability
p_unbinding=0.0
# grab radius - [A]
r_g=5.0E+2

# simulator parameters

# loop update frequency - [# timesteps]
freq_loop=10000

# topoisomerase relaxation (soft DNA-DNA pairs) frequency - [# timesteps]
freq_topo=50000
# topoisomerase relaxation (soft DNA-DNA pairs) interval - [# timesteps]
dNt_topo=50000
```
**Simulating with loops and topoisomerase:**\
In your second terminal, go to `/examples/simulating_chromosome_with_loops_and_topo/` and open 'simulating_chromosome_with_loops_and_topo.inp'. Make sure that "terminate" is commented, and "switch_skip_runs" is set to "F". 

In line 52, rename the LAMMPS simulation output to ",example_loops".

*Optional: You may consider increase the number of loops, for example by changing line 84 to`simulator_run_loops:30,100000,10000,25000,noappend,first`.

Save your changes and run the directives file. 

**Visualizing LAMMPS trajectories with VMD:**\
You will now copy over the .lammpstrj files from Delta to your local machine in order to visualize them in vmd.
```bash
scp $USERNAME@login.delta.ncsa.illinois.edu:/projects/bcuj/$USERNAME/btree_chromo_workspace/examples/simulating_chromosome_with_loops_and_topo/*.lammpstrj .
```
 In the Tk Console, do `source load_example_loops.tcl`. Anchors are colored black, while hinges are colored white.
 
## Links
Link to Software [bTreeChromo Github](https://github.com/brg4/btree_chromo)
Link to Pubilication  [Gilbert et al. Frontiers in Cell & Dev. Bio., 2023](https://www.frontiersin.org/articles/10.3389/fcell.2023.1214962/full)
