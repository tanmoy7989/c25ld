#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N convert
date

python convert_traj.py ~/c25ld/data/lammpstraj/methane_wca_transferability_runs/methane10_wca.lammpstrj.gz /home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane_wca_transferability_runs/methane10_wca_mapped.lammpstrj.gz 

python convert_traj.py ~/c25ld/data/lammpstraj/methane_wca_transferability_runs/methane20_wca.lammpstrj.gz /home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane_wca_transferability_runs/methane20_wca_mapped.lammpstrj.gz 

python convert_traj.py ~/c25ld/data/lammpstraj/methane_waterbox/methane30_wca.lammpstrj.gz /home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane_waterbox/methane30_wca_mapped.lammpstrj.gz 
