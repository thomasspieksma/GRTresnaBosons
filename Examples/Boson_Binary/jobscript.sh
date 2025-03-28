#!/bin/sh                                                  
#SBATCH -J Scalar_BH_2                    # job name                  
#SBATCH -t 1:00:00                      # walltime (dd:hh:mm:ss)    
#SBATCH -p astro2_devel

#SBATCH --nodes 8
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=4
# Output files
#SBATCH -o std_output.txt
#SBATCH -e std_error.txt
#SBATCH -D ./
#SBATCH --exclusive

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# Run the program
srun ./Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex params.txt

