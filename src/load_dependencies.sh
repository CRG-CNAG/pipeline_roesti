#!/bin/bash

# export LC_ALL="en_US.UTF-8"

# # DRMAA library for using Ruffus python pipeline on the SGE cluster
export DRMAA_LIBRARY_PATH=/usr/share/univage/lib/lx-amd64/libdrmaa.so

source /users/lserrano/mweber/Research_cloud/Python_mwTools/load_python_path_custom_modules.sh

# On simba
# source /home/mweber/.local/anaconda3/bin/activate roesti
# On cluster
source /users/lserrano/mweber/.local/miniconda3/bin/activate roesti
