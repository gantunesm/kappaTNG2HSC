#! /bin/bash

#SBATCH --nodes=1
#SBATCH --tasks-per-node=68
###32
#SBATCH --cpus-per-task=1
#SBATCH -C knl
##haswell
#SBATCH -t 04:00:00
#SBATCH -J dark_stats
#SBATCH -o dark_stats.o%j
#SBATCH -e dark_stats.e%j
#SBATCH --qos=regular
##regular
#SBATCH -A mp107
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gabriela.antunesm@gmail.com

#source /global/common/software/m3035/conda-activate.sh 3.7


python /global/homes/g/gmarques/hsc_ng/scripts/make_noisytng_dark_stats.py
