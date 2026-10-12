#!/bin/sh
# TestU01 batteries (SmallCrush, Crush, BigCrush) for V-MT19937-64 and V-SFMT19937,
# run in parallel on a remote Linux machine (one process per battery and configuration).
#
# Usage: ./run_testu01_remote.sh <start|status|collect> <ip> <user> <ssh-key>
#
#   ip       address or host name of the remote machine
#   user     user name on the remote machine
#   ssh-key  path of the private key used to connect
#
#   start    sync the repo, install TestU01 if missing, build, launch all the jobs in background
#   status   show which jobs are still running and the summary of the completed ones
#   collect  copy the logs into logs/testu01/
#
# The jobs keep running after the ssh session is closed. There are 21 jobs, of which 14
# (Crush and BigCrush) are long: a machine with 16 physical cores runs them all at once.
# TestU01 tests the output sequence, which does not depend on the hardware register width,
# so the generators are built with NBITS=128 and wider registers are emulated.
#
# Environment: NBITS (default 128),
#              TESTU01_ZIP (local copy of the TestU01 1.2.3 sources, default ../testu01/TestU01.zip)

set -e

if [ $# -ne 4 ]; then
    echo "Usage: $0 <start|status|collect> <ip> <user> <ssh-key>" >&2
    exit 1
fi

CMD="$1"
REMOTE="$3@$2"
SSH_KEY="$4"

REMOTE_DIR="vmt_testu01"
SSH_OPTS="-i $SSH_KEY -o StrictHostKeyChecking=no"
NBITS="${NBITS:-128}"
TESTU01_ZIP="${TESTU01_ZIP:-../testu01/TestU01.zip}"
LOG_DIR="logs/testu01"

# generator:nbits. A single state (mt64:64, sfmt:128) is the base generator, used as reference.
CONFIGS="mt64:64 mt64:128 mt64:256 mt64:512 sfmt:128 sfmt:256 sfmt:512"

case "$CMD" in
start)
    echo "=== Syncing repo to $REMOTE:$REMOTE_DIR ==="
    rsync -az --delete \
        --exclude='.git' \
        --exclude='.vs' \
        --exclude='bin-*' \
        --exclude='logs/' \
        --exclude='tu01-logs/' \
        --exclude='dat/matrix/' \
        -e "ssh $SSH_OPTS" \
        . "$REMOTE:$REMOTE_DIR"

    # the TestU01 sources are uploaded from the local copy, if there is one
    if [ -f "$TESTU01_ZIP" ]; then
        echo "=== Uploading $TESTU01_ZIP ==="
        ssh $SSH_OPTS "$REMOTE" "mkdir -p testu01"
        scp $SSH_OPTS "$TESTU01_ZIP" "$REMOTE:testu01/TestU01.zip"
    fi

    echo "=== Installing TestU01 (if missing), building and launching on $REMOTE ==="
    ssh $SSH_OPTS "$REMOTE" "REMOTE_DIR='$REMOTE_DIR' NBITS='$NBITS' CONFIGS='$CONFIGS' bash -s" <<'EOS'
set -e

# build tools
missing=""
for tool in g++ make unzip curl; do
    command -v $tool > /dev/null 2>&1 || missing="$missing $tool"
done
if [ -n "$missing" ]; then
    echo "installing missing tools:$missing"
    if command -v apt-get > /dev/null 2>&1; then
        sudo apt-get update -y -q
        sudo apt-get install -y -q build-essential unzip curl
    else
        sudo dnf install -y gcc gcc-c++ make unzip
    fi
fi
g++ --version | head -n1

# TestU01 1.2.3: uploaded copy, or else from the web site of its authors (Universite de Montreal)
T="$HOME/testu01"
if [ ! -f "$T/install/lib/libtestu01.a" ]; then
    echo "building TestU01 in $T (a few minutes)"
    mkdir -p "$T"
    cd "$T"
    if [ ! -f TestU01.zip ]; then
        curl -fsSLO http://simul.iro.umontreal.ca/testu01/TestU01.zip             || { echo "cannot download TestU01.zip: set TESTU01_ZIP to a local copy"; exit 1; }
    fi
    unzip -q -o TestU01.zip
    cd TestU01-1.2.3
    if ! (./configure --prefix="$T/install" --disable-shared && make -j"$(nproc)" && make install) > "$T/build.log" 2>&1; then
        # recent compilers reject some of the old C constructs used by TestU01
        echo "default build failed, retrying with -std=gnu89"
        make distclean > /dev/null 2>&1 || true
        (./configure --prefix="$T/install" --disable-shared CFLAGS="-O2 -std=gnu89" && make -j"$(nproc)" && make install) > "$T/build.log" 2>&1 \
            || { tail -n 30 "$T/build.log"; echo "TestU01 build failed, see $T/build.log"; exit 1; }
    fi
fi

# the generators
cd "$HOME/$REMOTE_DIR"
if pgrep -x testu01 > /dev/null 2>&1; then
    echo "testu01 jobs are already running: not starting again"
    exit 1
fi
make -k -j"$(nproc)" NBITS="$NBITS" TESTU01_DIR="$T/install" > build.log 2>&1 || true
exe=$(ls bin-"$NBITS"-*-g++/testu01 bin-"$NBITS"-*-g++/testu01.exe 2> /dev/null | head -n1)
if [ -z "$exe" ]; then
    tail -n 30 build.log
    echo "build of testu01 failed, see $HOME/$REMOTE_DIR/build.log"
    exit 1
fi
cp "$exe" ./testu01
echo "built $exe"

# one background process per battery and configuration
mkdir -p tu01-logs
n=0
for cfg in $CONFIGS; do
    gen=${cfg%%:*}
    bits=${cfg##*:}
    m=0
    for battery in SmallCrush Crush BigCrush; do
        nohup ./testu01 -g="$gen" -b="$bits" -m="$m" > "tu01-logs/$battery-$gen-$bits.log" 2>&1 < /dev/null &
        m=$((m + 1))
        n=$((n + 1))
    done
done
echo "launched $n jobs on $(nproc) cpus"
EOS
    echo "=== Started. Check progress with: $0 status $2 $3 $4 ==="
    ;;

status)
    ssh $SSH_OPTS "$REMOTE" "cd $REMOTE_DIR && bash -s" <<'EOS'
echo "running jobs: $(pgrep -x -c testu01)"
for f in tu01-logs/*.log; do
    if grep -q 'Summary results' "$f"; then
        if grep -q 'All tests were passed' "$f"; then
            echo "DONE     $f: all tests were passed"
        else
            echo "DONE     $f: p-values outside [0.001, 0.9990]:"
            sed -n '/p-value$/,/All other tests were passed/p' "$f" | grep -E '^ *[0-9]+ ' | sed 's/^/             /'
        fi
    else
        echo "RUNNING  $f ($(grep -c 'CPU time used' "$f") tests completed)"
    fi
done
EOS
    ;;

collect)
    mkdir -p "$LOG_DIR"
    scp $SSH_OPTS "$REMOTE:$REMOTE_DIR/tu01-logs/*.log" "$LOG_DIR/"
    echo "=== Logs saved in $LOG_DIR/ ==="
    ;;

*)
    echo "Usage: $0 <start|status|collect> <ip> <user> <ssh-key>" >&2
    exit 1
    ;;
esac
