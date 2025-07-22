#!/bin/bash

# Loop through all text files in the snakeATAC output directory
for file in output/counts/*.summary; do
    #echo "Processing file: $file"
    
    # Extract "Assigned" value
    assigned=$(awk '$1 == "Assigned" {print $2}' "$file")
    
    # Calculate the sum of all counts (excluding the first line with "Status")
    total=$(awk 'NR > 1 {sum += $2} END {print sum}' "$file")
    
    # Calculate the ratio (Assigned / Total)
    if [[ $total -ne 0 ]]; then
        ratio=$(echo "scale=6; $assigned / $total" | bc)
        #echo "Ratio (Assigned/Total): $ratio"
	echo $ratio
    else
        echo "Total is 0, cannot compute ratio"
    fi
    
    echo
done
