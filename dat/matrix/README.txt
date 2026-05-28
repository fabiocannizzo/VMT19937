Transition Matrices
===================

The jump-ahead transition matrices are large and are stored in a separate branch (`jump-matrix`) to keep the main repository lean.

To retrieve and extract the matrices, run the following commands:

1. Fetch the jump-matrix branch:
   git fetch origin jump-matrix:jump-matrix

2. Restore the matrix files into your working directory:
   git restore --source=jump-matrix -- "dat/matrix/**/*.bits"
   git restore --source=jump-matrix -- "dat/matrix/**/*.7z"

3. Extract the .7z archives and delete them:
   # Using find (Linux/Cygwin/Git Bash):
   find dat/matrix -name "*.7z" -exec sh -c '7z e "$1" -o"$(dirname "$1")" -y && rm "$1"' _ {} \;

   # Or using PowerShell:
   Get-ChildItem dat/matrix -Filter *.7z -Recurse | ForEach-Object { 7z e $_.FullName "-o$($_.DirectoryName)" -y; Remove-Item $_.FullName }
