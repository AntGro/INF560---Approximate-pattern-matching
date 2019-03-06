#!/bin/bash -l

cores=(1 2 4)
nodes=(1 2)

#A MODIFIER
#FAIRE UNE LISTE DE PATTERNS A PARCOURIR
#tester databases

for n in ${cores[@]}; do
    for N in ${nodes[@]}; do
        salloc -n $n -N $N ./apm 0 dna_database "A" "A" "A" "A"
    done

done





