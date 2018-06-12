#!/bin/sh
rm results.out

#GENERATE input.dat
#N_INPUT=4
#for i in `seq 1 $N_INPUT`
#do
#echo $i
#grep -v "#" all_input.dat | head -$i  | tail -1 > input.dat
#EXEC CG code and collect the results
../../../bin/main.x > results.out

#echo "`cat input.dat`" " = " "`cat results.dat`" >> all_results.dat

#done


