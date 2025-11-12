#!/bin/bash
#SBATCH --account=harvey
#SBATCH -c 20 				# REQUEST ONE CORE
#SBATCH -t 1-01:00			# runtime
#SBATCH -p gpu_quad,gpu			# partition
#SBATCH --gres=gpu:rtx8000:1,vram:48G
#SBATCH --qos=gpuquad_qos
#SBATCH --mem=48G			# memory
#SBATCH -o KM49-20251021_%x_%j.out		# outputs
#SBATCH -e KM49-20251021_%x_%j.err		# errors

module load gcc/14.2.0
module load conda/miniforge3/24.11.3-0
module load python/3.13.1
module load cuda/12.8

conda activate kilosort

cd /home/biw842/SpikeSorting/

python sort_spikes_NPX_KS4_KMversion.py -m KM49-50 -d 251021_g0

