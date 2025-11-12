#!/bin/bash
#SBATCH --account=harvey
#SBATCH -c 20 				# REQUEST ONE CORE
#SBATCH -t 1-00:00			# runtime
#SBATCH -p gpu_quad,gpu			# partition
#SBATCH --gres=gpu:rtx8000:1,vram:48G
#SBATCH --qos=gpuquad_qos
#SBATCH --mem=48G			# memory
#SBATCH -o 251003093810_%x_%j.out		# outputs
#SBATCH -e 251003093810_%x_%j.err		# errors

module load gcc/14.2.0
module load conda/miniforge3/24.11.3-0
module load python/3.13.1
module load cuda/12.8

conda activate DEEPLABCUT

cd /home/biw842/DeepLabCut/

python DLC_analyze.py -d 251003093810

