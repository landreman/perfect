#!/bin/bash -l

# Standard output and error:
#SBATCH -o ./perfectJob.out.%j
#SBATCH -e ./perfectJob.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J perfect

# Queue (Partition):
#SBATCH --partition=express

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3

# Wall clock limit:
#SBATCH --time=00:10:00

module load hdf5-mpi
module load petsc-real 

export PATH=${PATH}:${HDF5_HOME}/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HDF5_HOME}/lib

srun ../../perfect -ksp_view -mat_mumps_icntl_4 2
