#!/bin/sh                                                  
#SBATCH -J Scalar_BH_2                    # job name                  
#SBATCH -t 6:00:00                      # walltime (dd:hh:mm:ss)    
#SBATCH -p astro2_short

#SBATCH --nodes 8
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=2
# Output files
#SBATCH -o std_output.txt
#SBATCH -e std_error.txt
#SBATCH -D ./

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# Run the program
srun ./Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPT.MPI.ex params.txt
