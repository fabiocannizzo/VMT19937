#!/bin/sh

filename=logs/perf/$(cat /proc/cpuinfo | grep "model name" | head -n1 | sed -e 's/model name.*[:] //g' -e 's/[(]R[)]//g' -e 's/[ ]/_/g')

rm -f "$filename"_*

for nbits in 128 256 512; do
    builddir="build-${nbits}"
    perf_exe="${builddir}/perf.exe"
    [ -f "$perf_exe" ] || perf_exe="${builddir}/perf"   # Linux has no .exe suffix
    if [ -f "$perf_exe" ]; then
        f="${filename}_${nbits}"
        echo "RESULTS $nbits" >> "$f"
        "$perf_exe" -n 10 --dir dat >> "$f"
    fi
done
