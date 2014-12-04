#/bin/bash
rm -fr output.dat time.dat
mpirun -np 8 ./test
sharpness=`cat output.dat|awk '{print $2}'|tail -n 5|tr '\n' ' '`
sharpness_baseline=`cat output-baseline.dat|awk '{print $2}'|tail -n 5|tr '\n' ' '`
echo "Sharpness:          $sharpness"
echo "Sharpness baseline: $sharpness_baseline"
time=`cat time.dat|awk '{print $2}'|tail -n 5|tr '\n' ' '`
time_baseline=`cat time-baseline.dat|awk '{print $2}'|tail -n 5|tr '\n' ' '`
echo "Time:          $time"
echo "Time baseline: $time_baseline"
