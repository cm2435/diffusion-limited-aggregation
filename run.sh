num_runs=32
#echo "parallel ./run ::: {1..$num_runs}"
echo "num runs : " $num_runs

seq 1 $num_runs | parallel --verbose -j24 ./run {} "force_vector_jump"