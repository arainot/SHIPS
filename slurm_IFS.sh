#!/bin/sh
#SBATCH --job-name=IFS11
#SBATCH --account=ivsusers
# SBATCH --time 1440
#SBATCH --output=/home/alan/Documents/Thesis/SPHERE/spectra/HD93403/mcmc11/stdout_ifs.log
#SBATCH --error=/home/alan/Documents/Thesis/SPHERE/spectra/HD93403/mcmc11/stderr_ifs.log
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
# SBATCH --mem=6000
#SBATCH --partition=normal
#SBATCH --qos=normal

declare -i k=25
declare -i kmax=39

while [ $k -lt $kmax ]; do
  srun echo "Wavelength: " $k
  srun python ships_ifs_MCMC_auto.py $k &
  k=k+1
done

srun echo "Done"
