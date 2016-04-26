#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N geom_analyse

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/lammpstraj/c25/unconstrained_lj/c25_unbiased_mapped.lammpstrj.gz "(lj,AA)" 25 1 1700 geom_prop/c25_lj_longtraj_7.8
python geom_analyse.py ./modtraj/c25_lj_longtraj_7.8/c25_lj_SP.lammpstrj.gz "(lj,SP)" 25 1 1700 geom_prop/c25_lj_longtraj_7.8
python geom_analyse.py ./modtraj/c25_lj_longtraj_7.8/c25_lj_SPLD.lammpstrj.gz "(lj,SPLD)" 25 1 1700 geom_prop/c25_lj_longtraj_7.8
python geom_analyse.py ./modtraj/c25_lj_longtraj_7.8/c25_lj_LD.lammpstrj.gz "(lj,LD)" 25 1 1700 geom_prop/c25_lj_longtraj_7.8

