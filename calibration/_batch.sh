#!/bin/tcsh
#SBATCH --job-name=testser
#SBATCH --partition=pAll
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task= 50
#SBATCH --mail-user=jalil.nourisa@hzg.de
#SBATCH --mail-type=END
#SBATCH --account=nourisa
#SBATCH --output=job.o%j
#SBATCH --error=job.e%j
module load applications/python/3.8
setenv OMP_NUM_THREADS 50
python3 calibration/calibration.py 50
