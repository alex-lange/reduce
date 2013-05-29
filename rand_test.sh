#!/bin/bash

num_tests=100
path=../data/domination/05-20_dom_rand/graphs

# for ((n=10; n<=100; n+=10 ))
for ((n=10; n<=40; n+=10 ))
do 
    for p in 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95
    do
	echo **************************
	echo $n $p ${num_tests} ${path}/n${n}_p${p}
	echo $n $p ${num_tests} | ./gen_gr -r ${path}/n${n}_p${p}
    done
done
