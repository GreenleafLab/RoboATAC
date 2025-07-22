#!/bin/bash

module load cuda/11.2.0 
module load cudnn/8.1.1.33

# https://stackoverflow.com/questions/34534513/calling-conda-source-activate-from-bash-script
eval "$(conda shell.bash hook)"
conda activate chrombp
module load system cairo
module load system pango

echo "Live"
"$@"