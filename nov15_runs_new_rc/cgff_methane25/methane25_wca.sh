#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N methane_rc
date
python cg.py wca 0 1 25 1700 /home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane25/nvt/methane_wca_unbiased_mapped.lammpstrj.gz methane25_wca 6.5 SPLD
python cg.py wca 0 1 25 1700 /home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane25/nvt/methane_wca_unbiased_mapped.lammpstrj.gz methane25_wca_highrc 8.0 
date
