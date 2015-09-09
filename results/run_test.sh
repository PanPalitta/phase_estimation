#/bin/bash
rm -fr output.dat time.dat
mpirun -np 8 ../src/phase_estimation ../test.cfg
sharpness=`cat test_output.dat|awk '{print $2}'|tail -n 5|tr '\n' ' '`
sharpness_baseline=`cat baseline-output.dat|awk '{print $2}'|tail -n 5|tr '\n' ' '`
echo "Sharpness:          $sharpness"
echo "Sharpness baseline: $sharpness_baseline"
time=`cat test_time.dat|awk '{print $2}'|tail -n 5|tr '\n' ' '`
time_baseline=`cat baseline-time.dat|awk '{print $2}'|tail -n 5|tr '\n' ' '`
echo "Time:          $time"
echo "Time baseline: $time_baseline"
rm test_*dat
