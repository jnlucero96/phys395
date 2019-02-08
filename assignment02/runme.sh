#!/bin/bash

fbuild () {
    gfortran -llapack -fdefault-real-8 -O3 -fimplicit-none -Wall $@
}

echo "Hello. Beginning production run of Assignment 02."
echo "Building the utilities..."
fbuild -c utilities.f90
echo "Utilities built."
echo

echo "Now compiling problem 1..."
fbuild -c problem1.f90
fbuild -o problem1 utilities.o problem1.o

echo "Compiling problem 1 done. Running problem 1..."
./problem1 < DATA_FILES/DATA.dat > problem1.dat

echo "Running problem 1 finished. Answer acquired is:"
cat problem1.dat

echo

echo "Now compiling problem 2 ..."
fbuild -c problem2.f90 && fbuild -o problem2 utilities.o problem2.o

echo "Compiling problem 2 done. Running problem 2..."
./problem2 < DATA_FILES/DATA.dat > problem2.dat

echo "Running problem 2 finished. Answer acquired is:"
cat problem2.dat

echo "Making a pretty plot..."

python plot.py problem2

echo

echo "Now compiling problem 3..."
fbuild -c problem3.f90 && fbuild -o problem3 utilities.o problem3.o

echo "Compiling problem 3 done. Running problem 3..."
./problem3 < DATA_FILES/DATA.dat > problem3.dat

echo "Running problem 3 finished. Answer acquired is:"
cat problem3.dat

echo

echo "Now compiling problem 4..."
sed -ie "10s/3/7/" problem2.f90
sed -ie "9s/3/7/" problem3.f90
rm *.f90e

fbuild -c problem2.f90 && fbuild -o problem4_2 utilities.o problem2.o
fbuild -c problem3.f90 && fbuild -o problem4_3 utilities.o problem3.o

echo "Compiling problem 4 done. Running problem 4..."
./problem4_2 < DATA_FILES/DATA.dat > problem4.dat
./problem4_3 < DATA_FILES/DATA.dat >> problem4.dat

echo "Running problem 4 finished. Answer acquired is:"
cat problem4.dat

echo "Making plots..."
python plot.py problem4

# restor the file to the original state
sed -ie "10s/7/3/" problem2.f90
sed -ie "9s/7/3/" problem3.f90
rm *.f90e

echo

echo "Now compiling problem 5..."
fbuild -c problem5.f90
fbuild -o problem5 utilities.o problem5.o

echo "Compiling problem 5 done. Running problem 5 for n = 10..."
./problem5 < DATA_FILES/DATA.dat

echo "Running problem 5 finished."

echo "Making plots..."
python plot.py problem5

echo

# we are finished with the executables. Delete all of them
echo "Done with executables. Removing them from EXISTENCE..."
rm problem? problem?_? *.o

# we are done with the data files. Delete all of them
echo "Done with the data files. Removing them from EXISTENCE..."
rm *.dat

echo "Finished executing. Have a great day!"