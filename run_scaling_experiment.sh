#!/bin/sh
# Controlled SIMD-width scaling experiment.
# Builds V-family generators at 128/256/512-bit on a single remote AVX-512 machine
# and collects throughput results for all three widths.
#
# Usage: ./run_scaling_experiment.sh <ip> <user> <ssh-key>
#
#   ip       address or host name of the remote machine
#   user     user name on the remote machine
#   ssh-key  path of the private key used to connect
#
# Environment: CXX (remote compiler, default g++), MIN_REPS (minimum repetitions per experiment, default 100)

set -e

if [ $# -ne 3 ]; then
    echo "Usage: $0 <ip> <user> <ssh-key>" >&2
    exit 1
fi

REMOTE="$2@$1"
SSH_KEY="$3"
REMOTE_DIR="/tmp/vmt_scaling"
SSH_OPTS="-i $SSH_KEY -o StrictHostKeyChecking=no"
LOG_DIR="logs/perf/scaling"
CXX="${CXX:-g++}"
MIN_REPS="${MIN_REPS:-100}"

echo "=== Syncing repo to $REMOTE:$REMOTE_DIR ==="
rsync -az --delete \
    --exclude='.git' \
    --exclude='bin-*' \
    --exclude='logs/' \
    -e "ssh $SSH_OPTS" \
    . "$REMOTE:$REMOTE_DIR"

mkdir -p "$LOG_DIR"

for nbits in 128 256 512; do
    echo ""
    echo "=== Building NBITS=$nbits on $REMOTE ==="
    ssh $SSH_OPTS "$REMOTE" "
        set -e
        cd $REMOTE_DIR
        make clean NBITS=$nbits CXX=$CXX
        make NBITS=$nbits CXX=$CXX perf
    "

    echo "=== Running perf NBITS=$nbits ==="
    ssh $SSH_OPTS "$REMOTE" "
        set -e
        cd $REMOTE_DIR
        perf_exe=\$(ls bin-${nbits}-*-$CXX/perf.exe bin-${nbits}-*-$CXX/perf 2>/dev/null | head -n1)
        \"\$perf_exe\" --dir dat --no-original --no-stl --no-mkl --no-xmt --min-reps=$MIN_REPS --alpha=0.005 > /tmp/scaling_${nbits}.log
    "

    logfile="$LOG_DIR/scaling_${nbits}.log"
    scp $SSH_OPTS "$REMOTE:/tmp/scaling_${nbits}.log" "$logfile"
    echo "  saved -> $logfile"
done

echo ""
echo "=== Done. Results in $LOG_DIR/ ==="
