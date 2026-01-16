#!/bin/sh                                                  
#SBATCH -J Scalar_BH_2                    # job name                  
#SBATCH -t 100:00:00                      # walltime (dd:hh:mm:ss)    
#SBATCH -p astro3_long

#SBATCH --nodes 6
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G
# Output files
#SBATCH -o std_output.txt
#SBATCH -e std_error.txt
#SBATCH -D ./
#SBATCH --exclusive

#export UCX_NET_DEVICES=all
#export UCX_TLS=rc,self,sm

# replace UCX settings with:
export I_MPI_FABRICS=shm:ofi
unset UCX_NET_DEVICES
unset UCX_TLS
unset FI_PROVIDER

#export I_MPI_FABRICS=ucx
#unset I_MPI_OFI_PROVIDER
#unset FI_PROVIDER


# Run the program
srun ./Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex params.txt

