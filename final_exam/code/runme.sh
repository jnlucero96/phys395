#!/bin/bash


chmod 744 *.py # give execute privilege to python files

# prototype compile shorthand
function fbuild () {
    gfortran -llapack -fdefault-real-8 -fimplicit-none -Ofast $@
}

echo "====== Problem 1 BEGIN ======"
echo "Building problem 1..."
fbuild -o p1 problem1.f90

# output to data file
echo "Running problem 1..."
./p1 > p1_data.dat

# read data file to screen
echo "Answers found:"
cat p1_data.dat

echo "Cleaning up..."

echo "Finished..."
echo "====== Problem 1 END ======"

echo ""

echo "====== Problem 2 BEGIN ======"
echo "Building problem 2..."
fbuild -o p2 problem2.f90

# output to data file
echo "Running problem 2..."
./p2 > p2_data.dat

# read data file to screen
echo "Answers found:"
cat p2_data.dat

echo "Cleaning up..."

echo "Finished..."
echo "====== Problem 2 END ======"

echo ""

echo "====== Problem 3 BEGIN ======"
echo "Building problem 3..."
fbuild -o p3 problem3.f90

# output to data file
echo "Running problem 3..."
./p3 > p3_data.dat

echo "Plotting results..."
./plot_q3_solutions.py

echo "Cleaning up..."

echo "Finished..."
echo "====== Problem 3 END ======"

echo ""

echo "====== Problem 4 BEGIN ======"
echo "Building problem 4..."
fbuild -o p4 problem4.f90

# output to data file
echo "Running problem 4..."
./p4 > p4_data.dat

echo "Plotting results..."
./plot_q4_solutions.py

echo "Cleaning up..."

echo "Finished..."
echo "====== Problem 4 END ======"

echo ""

echo "====== Problem 5 BEGIN ======"
echo "Compiling..."
fbuild -o p5 problem5.f90
echo "Running..."
./p5 > p5_output.dat
echo "First 10 Lowest Energy Eigenvalues of Quantum Oscillator are:"
cat fort.11
echo "Plotting the energy eigenvectors..."
gnuplot --persist << EOF
reset
unset key
set title "Eigenstates of the Quantum Oscillator"
set xrange [-1:1]
set terminal pdf
set output "p5_wavefunctions_figure.pdf"
plot for [ii=0:9] "p5_output.dat" index ii u 1:2 w l
EOF

echo "Plotting the probability distributions..."
gnuplot --persist << EOF
reset
unset key
set title "Probability Density of Quantum Oscillator"
set xrange [-1:1]
set terminal pdf
set output "p5_densities_figure.pdf"
plot for [ii=0:9] "p5_output.dat" index ii u 1:3 w l
EOF

echo "Cleaning up..."

echo "Finished..."
echo "====== Problem 5 END ======"

echo "Removing extraneous files FROM EXISTENCE..."
rm p? *.dat fort.*

echo "Done. Have a nice day!"