#!/usr/bin/env bash

# Name of your script
script="gpa_array_gradientshap_alvis_param.sh"

# Job array parameters
delta_fields=1
initial_field=24000
final_field=25000
num_cases=$(( (final_field - initial_field) / delta_fields ))

echo "Submitting job array with $num_cases tasks"

# Submit as a job array
sbatch --array=1-$num_cases%50 $script $initial_field $delta_fields
# --array=1-$num_cases: run from 1 to $num_cases
# %50 run at most 50 at the same time
# $script script to run
# $initial_field initial index
# $delta_fields number of fields
