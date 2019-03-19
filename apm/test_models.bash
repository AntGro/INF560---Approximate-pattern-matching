#!/bin/bash -l

cores=(1 2 4 10 20)
nodes=(1 2 5 10)
patterns=('A' 'ACGTTA' 'G' 'AAAC' 'T' 'CT' 'C' 'GTT' 'AGTC' 'CGGTA' 'AA' 'ACCGGT' 'CG' 'ACGGTAC' 'ACT' 'CT' 'ACGTAT' 'AAA' 'ACCCGTTGGAT' 'AAC')
#A MODIFIER
#FAIRE UNE LISTE DE PATTERNS A PARCOURIR
#tester databases

for n in ${cores[@]}; do
    for N in ${nodes[@]}; do
        if [N>n]; then
            for m in `seq 1 20`; do            
                tmp = ("${patterns[@]:0:m}")
                salloc -n $n -N $N mpirun ./apm 0 dna/chr1_136K_alt.fa "${tmp[@]}" 
            done
        fi
    done
done





