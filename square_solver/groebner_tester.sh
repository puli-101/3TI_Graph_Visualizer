#!/bin/bash

#### Computes statistics about the system of equations solver

# Clear previous output
OUTPT_GRBNR=output_groebner.csv
OUTPT_RES=raw_time_info.txt


# MODIFY: Parameters
N=5
M=5
K=5
Q=17

echo "echo1". $OUTPT_GRBNR
echo "echo2". $OUTPT_RES

>> $OUTPT_GRBNR
2>> $OUTPT_RES

# Loop 1000 times
for i in {1..10000}; do
    # Run the Python script with arguments and append the output to output.txt
    #set max execution priority & set PATH
    #display exec time and max
    #   #setup environment
    echo "Running test number " . $i
    sudo env PATH=$PATH nice -n -20 /usr/bin/time -v /home/david/mambaforge/envs/isl_sage/bin/python3.11 groebner_solver.py -n=$N -m=$M -k=$K -q=$Q --csv 2> $OUTPT_RES >> $OUTPT_GRBNR
    echo -n $i, $N, $M, $K, $Q, >> time_stats.csv
    python3 tester_aux.py $OUTPT_RES >> time_stats.csv
done

echo $'\a'
echo $'\a'
echo $'\a'
