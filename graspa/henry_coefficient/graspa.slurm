#!/bin/bash
#SBATCH -J graspa
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -p gpu
#SBATCH --gres=gpu:5 #分配显卡数量

echo START: `date`
module load graspa
# gRASPA > output.txt
python main_henry_coefficient.py
echo FINISH: `date`
