#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N cg_methane
date

python md.py wca 0 1 50 1700  /home/cask0/home/tsanyal/c25ld/data/lammpstraj/nov15_methane_wca_transferability_runs/methane50_wca_mapped.lammpstrj.gz methane50_wca_SP /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_methane25/methane25_wca_SP 6.5

python md.py wca 0 1 50 1700  /home/cask0/home/tsanyal/c25ld/data/lammpstraj/nov15_methane_wca_transferability_runs/methane50_wca_mapped.lammpstrj.gz methane50_wca_SPLD /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_methane25/methane25_wca_SPLD 6.5

python md.py wca 0 1 50 1700  /home/cask0/home/tsanyal/c25ld/data/lammpstraj/nov15_methane_wca_transferability_runs/methane50_wca_mapped.lammpstrj.gz methane50_wca_LD /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_methane25/methane25_wca_LD 6.5

python md.py wca 0 1 50 1700  /home/cask0/home/tsanyal/c25ld/data/lammpstraj/nov15_methane_wca_transferability_runs/methane50_wca_mapped.lammpstrj.gz methane50_wca_highrc_SP /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_methane25/methane25_wca_highrc_SP 6.5

python md.py wca 0 1 50 1700  /home/cask0/home/tsanyal/c25ld/data/lammpstraj/nov15_methane_wca_transferability_runs/methane50_wca_mapped.lammpstrj.gz methane50_wca_highrc_SPLD /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_methane25/methane25_wca_highrc_SPLD 6.5

python md.py wca 0 1 50 1700  /home/cask0/home/tsanyal/c25ld/data/lammpstraj/nov15_methane_wca_transferability_runs/methane50_wca_mapped.lammpstrj.gz methane50_wca_highrc_LD /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_methane25/methane25_wca_highrc_LD 6.5
