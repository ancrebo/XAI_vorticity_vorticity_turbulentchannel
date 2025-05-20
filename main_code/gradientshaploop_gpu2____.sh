# Name of your script
script="gpa_gradientshap_berzelius_param.sh"
script_cont="gpa_gradientshap_berzelius_param.sh"

# Number of times to run the script
yeswait=2
nowait=0
num_runs=1
delta_fields=1
initial_field=25800
final_field=28000
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
        elif [ $i -eq 1 ]; then
            prev_job_id=$(sbatch --dependency=after:$prev_job_id $script $fini $ffin $yeswait | awk '{print $4}')
            echo "Initial wait"
        else
            # Submit subsequent jobs with a dependency on the previous job
            prev_job_id=$(sbatch --dependency=afterany:$prev_job_id $script_cont $prev_job_id $fini $ffin | awk '{print $4}')
            echo "Continue"
        fi
        echo "Submitted run $i with job ID $prev_job_id"
    done
done
