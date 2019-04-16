#!/bin/bash

# Load dependencies for roesti pipeline on the CRG cluster
# See cluster web pages for more information http://www.linux.crg.es
# echo "Exporting paths..."

# export PATH="/users/lserrano/mweber/bin/roesti:$PATH"

# export LC_ALL="en_US.UTF-8"

# # DRMAA library for using Ruffus python pipeline on the SGE cluster
export DRMAA_LIBRARY_PATH=/usr/share/univage/lib/lx-amd64/libdrmaa.so

# # Qt 5.5
# # Is needed for ngs-bits tools (SeqPurge)
# export LD_LIBRARY_PATH=/users/lserrano/mweber/Software/Qt5.5.0/5.5/gcc_64/lib/:$LD_LIBRARY_PATH
# export PATH="/users/lserrano/mweber/Software/Qt5.5.0/5.5/gcc_64/bin/:$PATH"

# # ngs-bits (includes SeqPurge)
# export PATH="/users/lserrano/mweber/Software/SeqPurge/ngs-bits/bin/:$PATH"

# # bowtie2
# # Symbolic link to the binaries in /users/lserrano/mweber/bin/roesti
# #export PATH=/users/lserrano/mweber/Software/bowtie2-2.2.9/:$PATH

# # SAMtools version 1.3.1
# # Symbolic link to the binaries in /users/lserrano/mweber/bin/roesti
# #export PATH="/users/lserrano/mweber/Software/samtools-1.3.1":$PATH

# # zlib (needed for bedtools)
# export LD_LIBRARY_PATH="/users/lserrano/mweber/Software/zlib-1.2.8/lib:$LD_LIBRARY_PATH"

# # Bedtools 2.26
# # Symbolic link to the binaries in /users/lserrano/mweber/bin/roesti
# #export PATH="/users/lserrano/mweber/Software/bedtools2/bin":$PATH

# # Loading Python virtual env
# source /users/lserrano/mweber/Research_cloud/Python_mwTools/load_python_virtualenv_cluster.sh

source /users/lserrano/mweber/Research_cloud/Python_mwTools/load_python_path_custom_modules.sh

# On simba
# source /home/mweber/.local/anaconda3/bin/activate roesti
# On cluster
source /users/lserrano/mweber/.local/miniconda3/bin/activate roesti
