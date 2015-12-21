#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N geom_analyse


python geom_analyse.py ./modtraj/c25/c25_wca_SP.lammpstrj.gz "(wca,SP)" 25 1 1700 geom_prop/c25 
python geom_analyse.py ./modtraj/c25/c25_wca_SPLD.lammpstrj.gz "(wca,SPLD)" 25 1 1700 geom_prop/c25 
python geom_analyse.py ./modtraj/c25/c25_wca_LD.lammpstrj.gz "(wca,LD)" 25 1 1700 geom_prop/c25 

python geom_analyse.py ./modtraj/c10/c10_wca_SP.lammpstrj.gz "(wca,SP)" 10 1 1700 geom_prop/c10 
python geom_analyse.py ./modtraj/c10/c10_wca_SPLD.lammpstrj.gz "(wca,SPLD)" 10 1 1700 geom_prop/c10 
python geom_analyse.py ./modtraj/c10/c10_wca_LD.lammpstrj.gz "(wca,LD)" 10 1 1700 geom_prop/c10 

python geom_analyse.py ./modtraj/c20/c20_wca_SP.lammpstrj.gz "(wca,SP)" 20 1 1700 geom_prop/c20 
python geom_analyse.py ./modtraj/c20/c20_wca_SPLD.lammpstrj.gz "(wca,SPLD)" 20 1 1700 geom_prop/c20 
python geom_analyse.py ./modtraj/c20/c20_wca_LD.lammpstrj.gz "(wca,LD)" 20 1 1700 geom_prop/c20 

python geom_analyse.py ./modtraj/c30/c30_wca_SP.lammpstrj.gz "(wca,SP)" 30 1 1700 geom_prop/c30 
python geom_analyse.py ./modtraj/c30/c30_wca_SPLD.lammpstrj.gz "(wca,SPLD)" 30 1 1700 geom_prop/c30 
python geom_analyse.py ./modtraj/c30/c30_wca_LD.lammpstrj.gz "(wca,LD)" 30 1 1700 geom_prop/c30 

python geom_analyse.py ./modtraj/c40/c40_wca_SP.lammpstrj.gz "(wca,SP)" 40 1 1700 geom_prop/c40 
python geom_analyse.py ./modtraj/c40/c40_wca_SPLD.lammpstrj.gz "(wca,SPLD)" 40 1 1700 geom_prop/c40 
python geom_analyse.py ./modtraj/c40/c40_wca_LD.lammpstrj.gz "(wca,LD)" 40 1 1700 geom_prop/c40 

python geom_analyse.py ./modtraj/c50/c50_wca_SP.lammpstrj.gz "(wca,SP)" 50 1 1700 geom_prop/c50 
python geom_analyse.py ./modtraj/c50/c50_wca_SPLD.lammpstrj.gz "(wca,SPLD)" 50 1 1700 geom_prop/c50 
python geom_analyse.py ./modtraj/c50/c50_wca_LD.lammpstrj.gz "(wca,LD)" 50 1 1700 geom_prop/c50 


date
