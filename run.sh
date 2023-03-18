num_runs=256
#echo "parallel ./run ::: {1..$num_runs}"
echo "num runs : " $num_runs

seq 1 $num_runs | parallel --verbose -j24 ./run {} "force_vector_jump"