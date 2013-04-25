#!/bin/bash

num_tests=100
path=./rand_tests/b

# for ((n=10; n<=100; n+=10 ))
for ((n=10; n<=40; n+=10 ))
do 
    for p in 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95
    do
	echo **************************
	echo $n $p ${num_tests} ${path}/n${n}_p${p}
	echo $n $p ${num_tests} | ./reduce -r ${path}/n${n}_p${p}
    done
done
