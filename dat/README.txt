Some biary files are availble as part of the repo, but they are in a separate branch.
To obtain them, execite the commands (from the repo base folder):
git fetch origin jumpfiles:jumpfiles
git checkout jumpfiles dat/

The files with extenssion "7z" in this directory have been compressed with the command:
7za a -t7z -m0=lzma2 filename.7z filename.bits
They can be decompressed using the makefile, with the command:
make matrix

The original binary file, which have extenstion ".bits", were saved using the function MT19937Matrix::toBin, passing as argument an ofstream opened with the flag ios::binary.

Jump files can be regenerated from srcatch using the jump.exe executable.
For example, if the folder dat/sfmt is empty, and we want to generate some jump files for the SFMT generator:
jump -g sfmt -p ./dat/sfmt -f 50 -j 10 -s 100
  -g sfmt      : sfmt generator
  -p dat/sfmt  : save files to this folder
  -f 50        : save all files for 2^n jumps, were n is a multiple of 50
  -j 10        : use 10 threads
  -s 100       : stop after saving the first files for 2^n jumps, where n>100
Note the utility will read the largest ".bits" fiule already contained in the directory and will use it as starting point for the calculations.




