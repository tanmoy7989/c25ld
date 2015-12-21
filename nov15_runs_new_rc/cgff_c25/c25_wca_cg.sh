#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N cg_md
date
#restarted from after SPLD
python cg.py wca 1 1 25 1700 /home/cask0/home/tsanyal/c25ld/data/lammpstraj/feb15_runs/unconstrained_wca/c25_unbiased_mapped.lammpstrj.gz c25_wca 7.8 SPLD

