#!/bin/bash -l

cores=(1 2 4 8 16 32)
nodes=(1 2 3 4 5 6)

#A MODIFIER
#FAIRE UNE LISTE DE PATTERNS A PARCOURIR
#tester databases

for n in cores ; do
    for N in nodes ; do
        srun -n n -N N ./apm 0 dna_database "AAG" "AGG"
    done

done





