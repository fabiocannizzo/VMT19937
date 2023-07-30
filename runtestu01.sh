# export PATH=/workspace/repos/testu01/install/bin:$PATH
bin-128/testu01 -b 32 -m 0 > testu01-logs/0-32.log &
bin-128/testu01 -b 32 -m 1 >  testu01-logs/1-32.log &
bin-128/testu01 -b 32 -m 2 >  testu01-logs/2-32.log &
bin-128/testu01 -b 128 -m 0 >  testu01-logs/0-128.log &
bin-128/testu01 -b 128 -m 1 >  testu01-logs/1-128.log &
bin-128/testu01 -b 128 -m 2 >  testu01-logs/2-128.log &
bin-128/testu01 -b 256 -m 0 >  testu01-logs/0-256.log &
bin-128/testu01 -b 256 -m 1 >  testu01-logs/1-256.log &
bin-128/testu01 -b 256 -m 2 >  testu01-logs/2-256.log &
bin-128/testu01 -b 512 -m 0 >  testu01-logs/0-512.log &
bin-128/testu01 -b 512 -m 1 >  testu01-logs/1-512.log &
bin-128/testu01 -b 512 -m 2 >  testu01-logs/2-512.log &
