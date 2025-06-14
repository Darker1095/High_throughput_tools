#!/bin/bash
#SBATCH -J graspa
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p gpu
#SBATCH --gres=gpu:1 

echo START: `date`
module load graspa
# gRASPA > output.txt
python main_mix_gcmc.py
echo FINISH: `date`
