num_runs=5
#echo "parallel ./run ::: {1..$num_runs}"
echo $num_runs
parallel ./run ::: $(seq 1 $num_runs)