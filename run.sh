num_runs=256
#echo "parallel ./run ::: {1..$num_runs}"
echo $num_runs
parallel -N 2 "./run {1} {2}" ::: $(seq 1 $num_runs) force_vector
