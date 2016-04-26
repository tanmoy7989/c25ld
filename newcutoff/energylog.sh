#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N energylog
date

python energylog.py wca 1 1 25 1700 /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/modtraj/c25/c25_wca_SPLD.lammpstrj.gz c25_wca_SPLD /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_c25/c25_wca_SPLD 7.8

python energylog.py wca 1 1 10 1700 /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/modtraj/c10/c10_wca_SPLD.lammpstrj.gz c10_wca_SPLD /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_c25/c25_wca_SPLD 7.8

python energylog.py wca 1 1 20 1700 /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/modtraj/c20/c20_wca_SPLD.lammpstrj.gz c20_wca_SPLD /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_c25/c25_wca_SPLD 7.8

python energylog.py wca 1 1 30 1700 /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/modtraj/c30/c30_wca_SPLD.lammpstrj.gz c30_wca_SPLD /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_c25/c25_wca_SPLD 7.8

python energylog.py wca 1 1 40 1700 /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/modtraj/c40/c40_wca_SPLD.lammpstrj.gz c40_wca_SPLD /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_c25/c25_wca_SPLD 7.8

python energylog.py wca 1 1 50 1700 /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/modtraj/c50/c50_wca_SPLD.lammpstrj.gz c50_wca_SPLD /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_c25/c25_wca_SPLD 7.8
