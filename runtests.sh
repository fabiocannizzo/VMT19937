#!/bin/sh

filename="results.txt"

rm $filename

make clean

for i in 128 256 512; do
   make NBITS=$i -j4
   echo RESULTS $i >> $filename
   ./bin-$i/perf.exe -n 10 >> $filename
done
