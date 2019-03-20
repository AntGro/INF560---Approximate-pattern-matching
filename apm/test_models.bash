#!/bin/bash -l

cores=(1 2 4 10 20)
nodes=(1 2 5 10)
patterns=('A' 'ACGTTA' 'G' 'AAAC' 'T' 'CT' 'C' 'GTT' 'AGTC' 'CGGTA' 'AA' 'ACCGGT' 'CG' 'ACGGTAC' 'ACT' 'CT' 'ACGTAT' 'AAA' 'ACCCGTTGGAT' 'AAC')
databases=('dna/chr6_GL000254v2_alt.fa')


for n in ${cores[@]}; do
    for N in ${nodes[@]}; do
        if [ $N -le $n ]; then
            for database in ${databases[@]}; do
               for m in `seq 1 19`; do
                    tmp=("${patterns[@]:0:m}")
                    echo  salloc -n $n -N $N mpirun ./apm 0 $database "${tmp[@]}"
                    salloc -n $n -N $N mpirun ./apm 0 $database "${tmp[@]}"
                done
            done
        fi
    done
done





