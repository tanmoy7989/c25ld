#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N cg_methane
date

# lj runs only for 6.5
#python md.py lj 0 1 25 1700  /home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane25/nvt/methane_lj_unbiased_mapped.lammpstrj.gz methane25_lj_SP /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_methane25/methane25_lj_SP 6.5

#python md.py lj 0 1 25 1700  /home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane25/nvt/methane_lj_unbiased_mapped.lammpstrj.gz methane25_lj_SPLD /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_methane25/methane25_lj_SPLD 6.5

#python md.py lj 0 1 25 1700  /home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane25/nvt/methane_lj_unbiased_mapped.lammpstrj.gz methane25_lj_LD /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_methane25/methane25_lj_LD 6.5
