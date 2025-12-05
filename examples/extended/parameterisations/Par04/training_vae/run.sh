#!/bin/bash
## To be run with 3 arguments defining the run name and the study name for wandb, and gpu_id

source /cvmfs/sft.cern.ch/lcg/views/LCG_105_cuda/x86_64-el9-gcc11-opt/setup.sh

nvidia-smi

pip install numpy h5py matplotlib scipy scikit-learn wandb tf2onnx onnxruntime

## Provide your wandb api key here
export WANDB_API_KEY=<WANDB-API-KEY>

mkdir -p validation checkpoint conversion generation
python /afs/cern.ch/user/p/praikwar/public/par04/training/train.py --run-name $1 --study-name $2 --gpu-ids $3

