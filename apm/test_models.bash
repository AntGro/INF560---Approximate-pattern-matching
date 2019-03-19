#!/bin/bash -l

cores=(1 2 4 10 20)
nodes=(1 2 5 10)
patterns=('A' 'ACGTTA' 'G' 'AAAC' 'T' 'CT' 'C' 'GTT' 'AGTC' 'CGGTA' 'AA' 'ACCGGT' 'CG' 'ACGGTAC' 'ACT' 'CT' 'ACGTAT' 'AAA' 'ACCCGTTGGAT' 'AAC')
databases=('dna/line_chrY.fa' 'dna/small_chrY.fa')


for n in ${cores[@]}; do
    for N in ${nodes[@]}; do
        if [ $n>=$N ]; then
            for database in ${databases[@]}; do
               for m in `seq 1 20`; do
                    tmp=("${patterns[@]:0:m}")
                    salloc -n $n -N $N mpirun ./apm 0 $database "${tmp[@]}"
                done
            done
        fi
    done
done





