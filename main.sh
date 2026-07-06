#!/bin/sh
#
# Submit script for nba
#
#SBATCH --account=ieor_lam         # Replace ACCOUNT with your group account name
#SBATCH --job-name=dro_m        # The job name.
#SBATCH -c 4                      # The number of cpu cores to use
#SBATCH -t 4-23:59                 # Runtime in D-HH:MM

pip install gurobipy
 
#Command to execute Python program
python dro_moment.py
 
#End of script