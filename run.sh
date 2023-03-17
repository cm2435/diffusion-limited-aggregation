num_runs=4
#echo "parallel ./run ::: {1..$num_runs}"
echo "num runs : " $num_runs

parallel --verbose -N 2 "./run {1} {2}" ::: $(seq 1 $num_runs) vanilla
