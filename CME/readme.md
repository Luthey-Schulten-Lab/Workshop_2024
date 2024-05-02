In the first day, you will launch several jobs to delta GPU node.

## Login into Delta login node

    ```bash
    ssh $USERNAME@delta.ncsa.illinois.edu
    ```

    You need to type you password and do 2FA.

##  Go into your projects directory and copy tutorials
    ```bash
    cd /projects/bcuj/$USERNAME
    # copy the source code for the workshop
    cp /projects/bcuj/twu4/Workshop_2024 ./
    # copy the apptainer
    cp /projects/bcuj/twu4/apptainer ./

    ```

## Launch the Jupyter Notebook

launch a juputer notebook on a delta GPU node using *srun* and ssh into the GPU node remotely to do the tutorials.

+ First: submit a job to delta GPU node
    ```bash
    srun --account=bcuj-delta-gpu --partition=gpuA100x4 --time=01:00:00 --mem=32g --gpus-per-node=1 --tasks-per-node=1 --cpus-per-task=8 --nodes=1 apptainer exec --nv --containall --bind /projects/bcuj/$USERNAME/Workshop_2024/:/workspace /projects/bcuj/$USERNAME/apptainer/qcb_workshop_2024.sif jupyter-notebook /workspace/ --no-browser --port=8811 --ip=0.0.0.0 --allow-root
    ```
    Then you should wait for Delta to allocate the resources for you, when you see something like this, it means you are good to proceed:
    ```bash
    srun: job 3546627 queued and waiting for resources
    srun: job 3546627 has been allocated resources
    WARNING: could not mount /etc/localtime: not a directory
    [I 19:07:57.203 NotebookApp] Writing notebook server cookie secret to /u/twu4/.local/share/jupyter/runtime/notebook_cookie_secret
    [I 19:07:58.314 NotebookApp] [jupyter_nbextensions_configurator] enabled 0.6.3
    [I 19:07:58.316 NotebookApp] Serving notebooks from local directory: /workspace
    [I 19:07:58.316 NotebookApp] Jupyter Notebook 6.4.12 is running at:
    [I 19:07:58.316 NotebookApp] http://gpua021.delta.ncsa.illinois.edu:8811/?token=b2e7ca15cd9dc3a6893a1273e359c88869225bc29d66c80c
    [I 19:07:58.316 NotebookApp]  or http://127.0.0.1:8811/?token=b2e7ca15cd9dc3a6893a1273e359c88869225bc29d66c80c
    [I 19:07:58.316 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
    [C 19:07:58.329 NotebookApp]

        To access the notebook, open this file in a browser:
            file:///u/twu4/.local/share/jupyter/runtime/nbserver-13-open.html
        Or copy and paste one of these URLs:
            http://gpua021.delta.ncsa.illinois.edu:8811/?token=b2e7ca15cd9dc3a6893a1273e359c88869225bc29d66c80c
        or http://127.0.0.1:8811/?token=b2e7ca15cd9dc3a6893a1273e359c88869225bc29d66c80c
    ```
+ Second: ssh into the delta GPU node.
    Your `$DeltaNode` can be found from the information above in last two lines after `http://`. For example, my Delta Node would be `gpua021`, all you need to do is replace `$DeltaNode` with your node you see above and replace `$USERNAME` with your username, for me its `twu4`.
    ```bash
    ssh -l $USERNAME  -L 127.0.0.1:8811:$DeltaNode.delta.internal.ncsa.edu:8811 dt-login.delta.ncsa.illinois.edu
    ```
    The command you should type in your local laptop console should be like this:
    ```bash
    ssh -l twu4 -L 127.0.0.1:8811:gpua021.delta.internal.ncsa.edu:8811 dt-login.delta.ncsa.illinois.edu
    ```

+ Third: open Jupyter Notebook in your own browser with the linke provided above. It will be something like this from the info shown by Delta:
    ``` bash
    http://127.0.0.1:8811/?token=b2e7ca15cd9dc3a6893a1273e359c88869225bc29d66c80c
    ```

## Launch CME/ODE Whole Cell Model

You will submit a job to delta GPU node using `SLRUM` system. 

+ First: go to folder

    ``` bash
    cd /projects/bcuj/$USERNAME/LM/CME/WCM/programs
    ```

+ Second: change the bash file and submit

    ```bash
    sbatch mpirun.sh
    ```

    


