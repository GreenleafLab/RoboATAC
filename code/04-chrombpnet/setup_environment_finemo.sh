#!/bin/bash

module load cuda/11.2.0 
module load cudnn/8.1.1.33

# https://stackoverflow.com/questions/34534513/calling-conda-source-activate-from-bash-script
eval "$(conda shell.bash hook)"
#conda activate finemo
conda activate finemo-dev # this dev branch supports using different peak sets for motif discovery and hit calling
module load system cairo
module load system pango

echo "Live"
"$@"