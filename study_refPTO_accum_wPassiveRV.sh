#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simmo536@umn.edu
#SBATCH -p small
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

cd ~/2023-NSF_CPS-wavePoweredRO
module load matlab
matlab -nodisplay -r \
"iVar = ${SLURM_ARRAY_TASK_ID}; \
study_refPTO_accum_wPassiveRV"

# Commands to use
# sbatch --array=1-675 ~/2023-NSF_CPS-wavePoweredRO/study_refPTO_accum_wPassiveRV.sh
# dos2unix  study_refPTO_accum_wPassiveRV.sh

