#!/usr/bin/bash

cd /home/jchishol/Summer_Project/
source /home/jchishol/opt/miniconda3/etc/profile.d/conda.sh
conda activate py2k
python plotter.py
conda deactivate
