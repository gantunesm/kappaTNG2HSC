#!/bin/bash 
#PBS -N stats_dark
#PBS -o /home/gmarques/ng_hsc/scripts
#PBS -e /home/gmarques/ng_hsc/scripts
#PBS -l select=1:ncpus=52
#PBS -l walltime=12:00:00
#PBS -u gmarques
#PBS -M gabriela.antunesm@gmail.com
#PBS -m ae
#PBS -q mini

# activate conda environment

source /home/anaconda3/bin/activate
 
conda activate tf39_cpu
 
python /home/gmarques/ng_hsc/scripts/make_noisytng_dark_stats.py
 	
