#!/bin/bash
#SBATCH --job-name=WholeCellModel_4replicates
#SBATCH --account=bbsv-delta-gpu
#SBATCH --output=WCM_4replicates.out
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=00:30:00
#SBATCH --mem=64g
#SBATCH --partition=gpuA100x4


export HYDRA_BOOTSTRAP=fork

# Run the Python script
apptainer exec --nv --bind /projects/bcuj/$USER/ /projects/bcuj/$USER/LM/LM.sif mpirun -np 4 python /projects/bcuj/$USER/LM/CME/WholeCellModel/programs/WCM_CMEODE_Hook.py -st cme-ode -t 120 -rs 60 -hi 1 -f "/projects/bcuj/"$USER"/LM/CME/WholeCellModel/output_4replicates"


# The bash file to launch parallel CMEODE simulations
# Each CMEODE simulation is independent with each other, i.e. do not communicate with each other

# This bash file will launch -np processes, each process execute python3 WCM_CMEODE_Hook.py -st cme-ode ...
# Each process may use single or multiple CPU cores/threads


# Arguments
# for mpirun:

    # -np numbers of instances/processes, integer number from 1 to nmax

# for python:

    # -st simulation type, only support "cme-ode"

    # -t simulation time, integer numbers, in seconds

    # -rs restart interval, integer numbers, in seconds

    # -hi hook interval, integer numbers, in seconds
    
    # -f directory to store output trajectory .csv files and log .txt files, strings, created automatically

    # For the times, the former should be the integer multiples of the latter e.g. -t 120 -rs 60 -wi 2 -hi 1
