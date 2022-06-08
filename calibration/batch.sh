#!/bin/tcsh
#SBATCH --job-name=testser
#SBATCH --partition=all
#SBATCH --nodes=4
#SBATCH --time=01:00:00
#SBATCH --mail-user=jalil.nourisa@hzg.de
#SBATCH --mail-type=ALL
#SBATCH --output=job.o%j
#SBATCH --error=job.e%j
unset LD_PRELOAD
# source /etc/profile.d/modules.sh
module purge

# module load applications/python/3.6
# setenv OMP_NUM_THREADS 20
python3 calibration/calibration.py 300
