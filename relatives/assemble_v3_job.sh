#!/bin/bash
#$ -cwd
#$ -j y
#$ -o logs/rel_assemble_v3.log
#$ -l h_data=4G,h_rt=4:00:00
#$ -pe shared 4
cd /u/project/pellegrini/gkislik/dogs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
./envs/genomics/bin/python3 relatives/assemble_v3.py
