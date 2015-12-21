#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N cg_methane
date

#wca runs for both cutoffs

python md.py wca 0 1 25 1700  /home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane25/nvt/methane_wca_unbiased_mapped.lammpstrj.gz methane25_wca_SP /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_methane25/methane25_wca_SP 6.5

python md.py wca 0 1 25 1700  /home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane25/nvt/methane_wca_unbiased_mapped.lammpstrj.gz methane25_wca_SPLD /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_methane25/methane25_wca_SPLD 6.5

python md.py wca 0 1 25 1700  /home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane25/nvt/methane_wca_unbiased_mapped.lammpstrj.gz methane25_wca_LD /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_methane25/methane25_wca_LD 6.5

python md.py wca 0 1 25 1700  /home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane25/nvt/methane_wca_unbiased_mapped.lammpstrj.gz methane25_wca_highrc_SP /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_methane25/methane25_wca_highrc_SP 8.0

python md.py wca 0 1 25 1700  /home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane25/nvt/methane_wca_unbiased_mapped.lammpstrj.gz methane25_wca_highrc_SPLD /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_methane25/methane25_wca_highrc_SPLD 8.0

python md.py wca 0 1 25 1700  /home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane25/nvt/methane_wca_unbiased_mapped.lammpstrj.gz methane25_wca_highrc_LD /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_methane25/methane25_wca_highrc_LD 8.0
