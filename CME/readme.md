## Code you will copy and run during the CME and RDME Tutorials.

## Login into Delta login node

```bash
ssh $USERNAME@login.delta.ncsa.illinois.edu
```

***Replace*** the $USERNAME with your username. You need to type you password and do 2FA.

##  Go into your projects directory and copy tutorials
```bash
# go to your directory
cd /projects/bcuj/$USER
```

```bash
# copy the source code for the workshop
cp -r /projects/bcuj/sharefile/Workshop_2024/LM ./
```

## Jupyter Notebook

launch a juputer notebook on a delta GPU node using *srun* and ssh into the GPU node remotely to do the tutorials.

+ First: submit a job to delta GPU node

    ***Replace*** $Port with a non-trivial number (Don't use 8888) to avoid using the same port as others.
    ```bash
    srun --account=bcuj-delta-gpu --partition=gpuA100x4 --time=08:00:00 --mem=64g --gpus-per-node=1 --tasks-per-node=1 --cpus-per-task=16 --nodes=1 apptainer exec --nv --containall --bind /projects/bcuj/$USER/:/workspace /projects/bcuj/$USER/LM/LM.sif jupyter-notebook /workspace/ --no-browser --port=$Port --ip=0.0.0.0 --allow-root
    ```
        ```bash
    srun --account=bcuj-delta-gpu --partition=gpuA100x4 --time=08:00:00 --mem=64g --gpus-per-node=1 --nodes=1 apptainer exec --nv --containall --bind /projects/bcuj/$USER/:/workspace /projects/bcuj/$USER/LM/LM.sif jupyter-notebook /workspace/ --no-browser --port=$Port --ip=0.0.0.0 --allow-root
    ```
    
    Then you should wait for Delta to allocate the resources for you, when you see something like this, it means you are good to proceed:
    ```bash
    srun: job 3546627 queued and waiting for resources
    srun: job 3546627 has been allocated resources
    WARNING: could not mount /etc/localtime: not a directory
    [I 19:07:57.203 NotebookApp] Writing notebook server cookie secret to /u/$USER/.local/share/jupyter/runtime/notebook_cookie_secret
    [I 19:07:58.314 NotebookApp] [jupyter_nbextensions_configurator] enabled 0.6.3
    [I 19:07:58.316 NotebookApp] Serving notebooks from local directory: /workspace
    [I 19:07:58.316 NotebookApp] Jupyter Notebook 6.4.12 is running at:
    [I 19:07:58.316 NotebookApp] http://gpua021.delta.ncsa.illinois.edu:8811/?token=b2e7ca15cd9dc3a6893a1273e359c88869225bc29d66c80c
    [I 19:07:58.316 NotebookApp]  or http://127.0.0.1:$Port/?token=b2e7ca15cd9dc3a6893a1273e359c88869225bc29d66c80c
    [I 19:07:58.316 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
    [C 19:07:58.329 NotebookApp]

        To access the notebook, open this file in a browser:
            file:///u/$USERNAME/.local/share/jupyter/runtime/nbserver-13-open.html
        Or copy and paste one of these URLs:
            http://`$DeltaNode`.delta.ncsa.illinois.edu:$Port/?token=b2e7ca15cd9dc3a6893a1273e359c88869225bc29d66c80c
        or http://127.0.0.1:$Port/?token=b2e7ca15cd9dc3a6893a1273e359c88869225bc29d66c80c
    ```

    The last two line contains the delta GPU node `$DeltaNode`.

+ Second: ssh into the delta GPU node.
    Open a second terminal.
  Your `$DeltaNode` can be found from the information above in last two lines after `http://`. ***replace*** `$DeltaNode` with your node you see above and ***replace*** `$USERNAME` with your username, for me its `twu4`. ***Replace*** `$Port` with the 4 digit number you used.
    
    ```bash
    ssh -l $USERNAME  -L 127.0.0.1:$Port:$DeltaNode.delta.internal.ncsa.edu:$Port dt-login.delta.ncsa.illinois.edu
    ```

    You need to type you password and do 2FA AGAIN.

+ Third: Copy the last line in the first terminal and paste to one browser to open Jupyter Notebook.

    ``` bash
    http://127.0.0.1:8811/?token=b2e7ca15cd9dc3a6893a1273e359c88869225bc29d66c80c
    ```

## CME/ODE Whole Cell simulation in Parallel

You will submit a job to run Whole Cell Model in parallel on delta GPU node.

In the given bash file, you will launch 2 minutes simulation of 4 replicates.

+ First: go to *programs* folder

    ``` bash
    cd /projects/bcuj/$USER/LM/CME/WholeCellModel/programs
    ```

+ Second: Just submit the bash file

    ```bash
    sbatch mpirun.sh
    ```
+ Third: Check your job
    Check the status of your job. *PD* means waiting to run, *R* running.

    ```bash
    squeue -u $USER
    ```
    go to *output_4replicates* folder 
    
    ``` bash
    cd /projects/bcuj/$USER/LM/CME/WholeCellModel/output_4replicates
    ```

