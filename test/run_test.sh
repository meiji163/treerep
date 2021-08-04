#! /bin/bash

declare -a tests=("immune" "randtree")
for t in ${tests[@]}; do
	$(g++ -std=c++11 -o ${t} ${t}.cc ../src/treerep.cc ../src/graph.cc) 
	echo "$(./${t})"
done
