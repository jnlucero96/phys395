#!/bin/bash

function fbuild () {
    gfortran -llapack -fdefault-real-8 -Ofast -fimplicit-none -lcfitsio $@
}

echo "Hello. Beginning production run of Assignment 04."

echo "Now compiling problem 1..."
fbuild -o p1 problem1.f90

echo "Compiling problem 1 done. Running problem 1..."
./p1 > problem1.dat

echo "Running problem 1 finished. Answer acquired is:"
cat problem1.dat

echo "Making plots..."
python plot.py q1_data.fit 1

echo

echo "Now compiling problem 2 ..."
fbuild -o p2 problem2.f90

echo "Compiling problem 2 done. Running problem 2..."
./p2 > problem2.dat

echo "Making plots..."
python plot.py problem2.dat 2

echo

echo "Now compiling problem 3..."
fbuild -o p3 problem3.f90

echo "Compiling problem 3 done. Running problem 3..."
./p3 > problem3.dat

echo "Running problem 3 finished. Answer acquired is:"
cat problem3.dat

echo "Making plots..."
python plot.py q3_data.fit 3

echo

# we are finished with the executables. Delete all of them
echo "Done with executables. Removing them from EXISTENCE..."
rm problem? *.o

# we are done with the data files. Delete all of them
echo "Done with the data files. Removing them from EXISTENCE..."
rm *.dat

echo "Finished executing..."

echo "If you want to animate the phase-space motion please run ./do_animate.gpl"
echo "Have a great day!"