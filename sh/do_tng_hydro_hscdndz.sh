#! /bin/bash

#SBATCH --nodes=1
#SBATCH --tasks-per-node=68
###32
#SBATCH --cpus-per-task=1
#SBATCH -C knl
##haswell
#SBATCH -t 03:00:00
#SBATCH -J make_hydro
#SBATCH -o make_hydro.o%j
#SBATCH -e make_hydro.e%j
#SBATCH --qos=regular
##regular
#SBATCH -A mp107
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gabriela.antunesm@gmail.com

#source /global/common/software/m3035/conda-activate.sh 3.7

#module load python3
python /global/homes/g/gmarques/hsc_ng/scripts/baryon_sims_hydro.py
