#!/bin/bash -l

cores=(1 2 4 10 20)
nodes=(1 2 5 10)

#A MODIFIER
#FAIRE UNE LISTE DE PATTERNS A PARCOURIR
#tester databases

for n in ${cores[@]}; do
    for N in ${nodes[@]}; do
        salloc -n $n -N $N mpirun ./apm 0  dna/chr1_136K_alt.fa AGTT AT AGTT AGTGTGT AGTCTC ACCCCCTG CCCCAT
    done
done





