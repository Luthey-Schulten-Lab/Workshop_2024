In the first day, you will launch several jobs to delta GPU node.

1. Login into Delta login node

ssh $USERNAME@delta.ncsa.illinois.edu

You need to type you password and do 2FA.

2. Go into your projects directory and copy tutorials

    cd /projects/bcuj/$USERNAME

    cp /projects/bcuj/enguang/LM ./


3. Launch Jupyter Notebook

   launch a juputer notebook on a delta GPU node using *srun* and ssh into the GPU node remotely to do the tutorials.

   First: submit a job to delta GPU node

   Second: ssh into the delta GPU node

   Third: open Jupyter Notebook in your own browser


4. Launch CME/ODE Whole Cell Model

    You will submit a job to delta GPU node using SLRUM system. 

    First: go to folder

    cd /projects/bcuj/$USERNAME/LM/CME/WCM/programs

    Second: change the bash file and submit

    sbatch mpirun.sh

    


