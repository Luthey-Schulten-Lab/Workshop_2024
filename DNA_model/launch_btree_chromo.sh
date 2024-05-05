#!/bin/bash

# set the location of the btree_chromo files
btree_chromo_files='/projects/bcuj/sharefile/Workshop_2024/DNA_model'

# set the main directory for the project/workshop and the user's subdirectory
project_dir='/projects/bcuj'
user_subdir=${project_dir}/${USER}

# use the workspace
workspace_dir=${user_subdir}/btree_chromo_workspace

# use srun to launch an interactive session
srun \
    --pty \
    --account=bbsv-delta-cpu \
    --partition=cpu \
    --time=00:30:00 \
    --mem=32g \
    --tasks-per-node=1 \
    --cpus-per-task=8 \
    --nodes=1 apptainer shell \
    --writable-tmpfs \
    --no-home \
    --containall \
    --bind ${workspace_dir}:/mnt\
    ${btree_chromo_files}/btree_chromo_2024workshop.sif
