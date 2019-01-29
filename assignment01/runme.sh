#!/bin/bash

echo "Hello. Beginning production run of Assignment 01."
echo

echo "Now compiling problem 1..."
gfortran problem1.f90 -fdefault-real-8 -O3 -fimplicit-none -Wall -l lapack -o problem1

echo "Compiling problem 1 done. Running problem 1..."
./problem1 > problem1.dat

echo "Running problem 1 finished. Answer acquired is:"
cat problem1.dat

echo

echo "Now compiling problem 2 for n = 10..."
gfortran problem2.f90 -fdefault-real-8 -O3 -fimplicit-none -Wall -l lapack -o problem2

echo "Compiling problem 2 done. Running problem 2 for n = 10..."
./problem2 > problem2_10n.dat

echo "Running problem 2 for n = 10 finished. Answer acquired is:"
cat problem2_10n.dat

echo

echo "Re-compiling problem 2 for n = 100..."
sed -ie "8s/n = 10/n = 100/" problem2.f90
rm problem2.f90e
gfortran problem2.f90 -fdefault-real-8 -O3 -fimplicit-none -Wall -l lapack -o problem2

echo "Compiling problem 2 for n = 100 done. Running problem 2 for n = 100..."
./problem2 > problem2_100n.dat
# restore the file to its original state after modifying it
sed -ie "8s/n = 100/n = 10/" problem2.f90
rm problem2.f90e

echo "Running problem 2 for n = 100 finished. Answer acquired is:"
cat problem2_100n.dat

echo

echo "Now compiling problem 3 for n = 10..."
gfortran problem3.f90 -fdefault-real-8 -O3 -fimplicit-none -Wall -l lapack -o problem3

echo "Compiling problem 3 done. Running problem 3 for n = 10..."
./problem3 > problem3_10n.dat

echo "Running problem 3 for n = 10 finished. Answer acquired is:"
cat problem3_10n.dat

echo

echo "Re-compiling problem 3 for n = 100..."
sed -ie "8s/n = 10/n = 100/" problem3.f90
rm problem3.f90e
gfortran problem3.f90 -fdefault-real-8 -O3 -fimplicit-none -Wall -l lapack -o problem3

echo "Compiling problem 3 for n = 100 done. Running problem 3 for n = 100..."
./problem3 > problem3_100n.dat
# restore the file to its original state after modifying it
sed -ie "8s/n = 100/n = 10/" problem3.f90
rm problem3.f90e

echo "Running problem 3 for n = 100 finished. Answer acquired is:"
cat problem3_100n.dat

echo

echo "Now compiling problem 4 for n = 10..."
gfortran problem4.f90 -fdefault-real-8 -O3 -fimplicit-none -Wall -l lapack -o problem4

echo "Compiling problem 4 done. Running problem 4 for n = 10..."
./problem4 > problem4_10n.dat

echo "Running problem 4 for n = 10 finished. Answer acquired is:"
cat problem4_10n.dat

echo

echo "Re-compiling problem 4 for n = 100..."
sed -ie "8s/n = 10/n = 100/" problem4.f90
rm problem4.f90e
gfortran problem4.f90 -fdefault-real-8 -O3 -fimplicit-none -Wall -l lapack -o problem4

echo "Compiling problem 4 for n = 100 done. Running problem 4 for n = 100..."
./problem4 > problem4_100n.dat
# restore the file to its original state after modifying it
sed -ie "8s/n = 100/n = 10/" problem4.f90
rm problem4.f90e

echo "Running problem 4 for n = 100 finished. Answer acquired is:"
cat problem4_100n.dat

echo

echo "Now compiling problem 5 for n = 10..."
gfortran problem5.f90 -fdefault-real-8 -O3 -fimplicit-none -Wall -l lapack -o problem5

echo "Compiling problem 5 done. Running problem 5 for n = 10..."
./problem5 > problem5_10n.dat

echo "Running problem 5 for n = 10 finished. Answer acquired is:"
cat problem5_10n.dat

echo

echo "Re-compiling problem 5 for n = 100..."
sed -ie "8s/n = 10/n = 100/" problem5.f90
rm problem5.f90e
gfortran problem5.f90 -fdefault-real-8 -O3 -fimplicit-none -Wall -l lapack -o problem5

echo "Compiling problem 5 for n = 100 done. Running problem 5 for n = 100..."
./problem5 > problem5_100n.dat
# restore the file to its original state after modifying it
sed -ie "8s/n = 100/n = 10/" problem5.f90
rm problem5.f90e

echo "Running problem 5 for n = 100 finished. Answer acquired is:"
cat problem5_100n.dat

echo

# we are finished with the executables. Delete all of them
echo "Done with executables. Removing them from EXISTENCE..."
rm problem? 

echo "Making pretty plots..."

python pretty_plots.py problem2
python pretty_plots.py problem3
python pretty_plots.py problem5

# we are done with the data files. Delete all of them
echo "Done with the data files. Removing them from EXISTENCE..."
rm *.dat

echo "Opening all plots..."
open *.pdf

echo "Finished executing. Have a great day!"