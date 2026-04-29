#!/bin/bash
# Benchmark automation script for Arch Linux (SSE4.2, AVX2, AVX-512)
# Assumes intel-oneapi-mkl is installed via pacman

# 1. Path Detection
MKL_BASE="/opt/intel/oneapi/mkl"
if [ -d "$MKL_BASE/latest" ]; then
    MKL_INC="-I$MKL_BASE/latest/include"
    MKL_LIB="-L$MKL_BASE/latest/lib"
else
    # Fallback to the first versioned directory found
    VERSION=$(ls $MKL_BASE | grep -E '^[0-9]' | head -n 1)
    if [ -n "$VERSION" ]; then
        MKL_INC="-I$MKL_BASE/$VERSION/include"
        MKL_LIB="-L$MKL_BASE/$VERSION/lib"
    else
        echo "Warning: MKL not found in $MKL_BASE. Reference tests may fail."
    fi
fi

LOG_DIR="logs/perf/linux-arch-results"
mkdir -p "$LOG_DIR"

echo "Starting VMT19937 Performance Suite..."
echo "MKL Include: $MKL_INC"

# 2. Benchmarking Loop
for NBITS in 128 256 512; do
    echo "------------------------------------------------"
    echo "Configuring ${NBITS}-bit SIMD..."
    BINDIR="bin-${NBITS}-g++"

    make clean > /dev/null
    make NBITS=$NBITS $BINDIR/perf.exe MKL_INC="$MKL_INC" MKL_LIB="$MKL_LIB"

    if [ -f "$BINDIR/perf.exe" ]; then
        echo "Running ${NBITS}-bit benchmark (slow mode)..."
        ./$BINDIR/perf.exe --slow > "$LOG_DIR/perf-${NBITS}.log" 2>&1
        echo "Results saved to $LOG_DIR/perf-${NBITS}.log"
    else
        echo "Error: Failed to build ${NBITS}-bit configuration."
    fi
done

# 3. Metadata Collection
echo "Collecting system specifications..."
{
  echo "=== CPU INFO ==="
  lscpu
  echo -e "\n=== COMPILER ==="
  g++ --version
  echo -e "\n=== MKL PACKAGE ==="
  pacman -Qi intel-oneapi-mkl 2>/dev/null || echo "MKL package info not available via pacman"
  echo -e "\n=== OS/KERNEL ==="
  uname -a
} > "$LOG_DIR/system_specs.txt"

echo "------------------------------------------------"
echo "All benchmarks complete."
echo "Final results located in: $LOG_DIR"
