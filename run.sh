num_runs=4
#echo "parallel ./run ::: {1..$num_runs}"
echo "num runs : " $num_runs

seq 1 $num_runs | parallel -j0 ./run {} "force_vector"