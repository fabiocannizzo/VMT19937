#!/bin/sh

filename=logs/perf/$(cat /proc/cpuinfo | grep "model name" | head -n1 | sed -e 's/model name.*[:] //g' -e 's/[(]R[)]//g' -e 's/[ ]/_/g')

rm -f "${filename}"_*

for nbits in 128 256 512; do
    perf_exe=$(ls bin-${nbits}-*-g++/perf.exe bin-${nbits}-*-g++/perf 2>/dev/null | head -n1)
    if [ -n "$perf_exe" ] && [ -f "$perf_exe" ]; then
        f="${filename}_${nbits}"
        echo "RESULTS $nbits" >> "$f"
        "$perf_exe" --dir dat >> "$f"
    fi
done
