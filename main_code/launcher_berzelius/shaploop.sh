# Name of your script
script="launch_gradientshap_tf2_11_h5.sh"

# Number of times to run the script
yeswait=5
nowait=0
num_runs=1
delta_fields=50
initial_field=20000
final_field=21000
num_cases=$(echo "scale=0; (($final_field - $initial_field) / $delta_fields) " | bc)
echo $num_cases

# Initialize the job ID variable
prev_job_id=""

# Loop to run multiple cases
for (( j=1; j<=num_cases; j++ )); do
    # Loop to submit the same script multiple times
    fini=$initial_field
    echo $fini
    echo "$j < $num_cases"
    if [ "$j" -lt "$num_cases" ]; then
        ffin=$(($initial_field+$delta_fields))
        initial_field=$ffin
    else
        ffin=$final_field
    fi
    echo $ffin
    for (( i=1; i<=num_runs; i++ )); do
        if [ -z "$prev_job_id" ]; then
            # Submit the first job without dependency
            prev_job_id=$(sbatch $script $fini $ffin $nowait | awk '{print $4}')
            echo "Initial"
        else
            prev_job_id=$(sbatch --dependency=after:$prev_job_id $script $fini $ffin $yeswait | awk '{print $4}')
            echo "Initial wait"
        fi
        echo "Submitted run $i with job ID $prev_job_id"
    done
done
