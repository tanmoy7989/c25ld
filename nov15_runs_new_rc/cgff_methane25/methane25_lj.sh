#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N methane_rc
date
python cg.py lj 0 1 25 1700 /home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane25/nvt/methane_lj_unbiased_mapped.lammpstrj.gz methane25_lj 6.5 
python cg.py lj 0 1 25 1700 /home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane25/nvt/methane_lj_unbiased_mapped.lammpstrj.gz methane25_lj_highrc 8.0 
date
