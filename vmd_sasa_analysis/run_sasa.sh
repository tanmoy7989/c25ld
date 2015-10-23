#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N vmd_sasa
date
python test_sasa.py 10 1700



