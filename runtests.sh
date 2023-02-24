#!/bin/sh

function test() {
   make clean
   make NBITS=$1 -j4
   echo RESULTS $1 >> $filename
   ./perf.exe -n 10 >> $filename
}

filename="results.txt"

rm $filename

test 128
test 256
test 512
