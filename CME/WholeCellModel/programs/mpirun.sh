mpirun -np 25 python3 WCM_CMEODE_Hook.py -st cme-ode -t 6300 -rs 60 -hi 1 -f '../output_10thApril_1' 

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
