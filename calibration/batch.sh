#!/bin/tcsh
#SBATCH --job-name=testser
#SBATCH --partition=pAll
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --mail-user=jalil.nourisa@hzg.de
#SBATCH --mail-type=END
#SBATCH --account=nourisa
#SBATCH --output=job.o%j
#SBATCH --error=job.e%j
./calibration.py 10
