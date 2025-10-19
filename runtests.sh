#!/bin/sh

#filename="results.txt"
filename=logs/perf/$(cat /proc/cpuinfo | grep "model name" | head -n1 | sed -e 's/model name.*[:] //g' -e 's/[(]R[)]//g' -e 's/[ ]/_/g')

rm $filename

make clean

for i in 128 256 512; do
    if [ -e ./bin-$i/perf.exe ]; then
        echo "File or directory exists."
      f=${filename}_$(NBITS}
      echo RESULTS $i >> ${f}
      ./bin-$i/perf.exe -n 10 >> ${f}
   fi
done
