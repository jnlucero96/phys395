
function fbuild ()
{
    gfortran -llapack -fdefault-real-8 -Ofast -fimplicit-none $@
}

echo "Starting Assignment #5..."

echo "====== Problem 1 BEGIN ======"
echo "Part A: even function"
echo "Compiling..."
fbuild -o p1_A problem1.f90
echo "Running..."
./p1_A > q1_output_odd.dat
echo "Has output:"
cat q1_output_odd.dat
echo "Plotting..."
gnuplot --persist << EOF
set title "Even Wavefunction (+ x-axis) Harmonic Oscillator"
set xrange [0:5]
set terminal pdf
set output "q1_even_function.pdf"
plot "fort.1" u 1:2 w l
EOF

echo ""

echo "Part B: odd function"
echo "Recompiling..."

sed -ie "15s/psi_init = 1.0/psi_init = 0.0/" problem1.f90
sed -ie "16s/psiPrime_init = 0.0/psiPrime_init = 1.0/" problem1.f90
sed -ie "17s/file_index_to_write = 1/file_index_to_write = 2/" problem1.f90
rm *.f90e

fbuild -o p1_B problem1.f90
echo "Running..."
./p1_B > q1_output_even.dat
echo "Has output:"
cat q1_output_even.dat
echo "Plotting..."
gnuplot --persist << EOF
set title "Odd Wavefunction (+ x-axis) Harmonic Oscillator"
set xrange [0:5]
set terminal pdf
set output "q1_odd_function.pdf"
plot "fort.2" u 1:2 w l
EOF

echo "Cleaning up..."

sed -ie "15s/psi_init = 0.0/psi_init = 1.0/" problem1.f90
sed -ie "16s/psiPrime_init = 1.0/psiPrime_init = 0.0/" problem1.f90
sed -ie "17s/file_index_to_write = 2/file_index_to_write = 1/" problem1.f90
rm *.f90e

echo "Finished..."
echo "====== Problem 1 END ======"

echo ""

echo "====== Problem 2 BEGIN ======"
echo "Compiling..."
fbuild -o p2 problem2.f90
echo "Running..."
./p2 > q2_output.dat
echo "First 10 Energy Eigenvalues of Harmonic Oscillator are:"
cat fort.3 
echo "Finished..."
echo "====== Problem 2 END ======"

echo ""

echo "====== Problem 3 BEGIN ======"
echo "Compiling..."

sed -ie "12s/option = 0/option = 1/" problem2.f90
rm *.f90e

fbuild -o p3 problem2.f90
echo "Running..."
./p3 > q3_output.dat
echo "First 10 Energy Eigenvalues of Quartic Oscillator are:"
cat fort.4

echo "Cleaning up..."

sed -ie "12s/option = 1/option = 0/" problem2.f90
rm *.f90e

echo "Finished..."
echo "====== Problem 3 END ======"

echo "====== Problem 4 BEGIN ======"
echo "Compiling..."
fbuild -o p4 problem4.f90
echo "Running..."
./p4 > q4_output.dat
echo "First 10 Energy Eigenvalues of Harmonic Oscillator are:"
cat fort.11
echo "Plotting the energy eigenvectors..."
gnuplot --persist << EOF
reset
unset key
set title "Eigenstates of the Harmonic Oscillator"
set xrange [-6:6]
set terminal pdf
set output "q4_wavefunctions_figure.pdf"
plot for [ii=0:9] "q4_output.dat" index ii u 1:2 w l
EOF

echo "Plotting the probability distributions..."
gnuplot --persist << EOF
reset
unset key
set title "Probability Density of the Harmonic Oscillator"
set xrange [-6:6]
set terminal pdf
set output "q4_densities_figure.pdf"
plot for [ii=0:9] "q4_output.dat" index ii u 1:3 w l
EOF

echo "Cleaning up..."

echo "Finished..."
echo "====== Problem 4 END ======"

echo "====== Problem 5 BEGIN ======"
echo "Compiling..."
sed -ie "8s/option = 0/option = 1/" problem4.f90
sed -ie "31s/file_out_index = 11/file_out_index = 12/" problem4.f90
rm *.f90e

fbuild -o p5 problem4.f90
echo "Running..."
./p5 > q5_output.dat
echo "First 10 Energy Eigenvalues of Quartic Oscillator are:"
cat fort.12
echo "Plotting the energy eigenvectors..."
gnuplot --persist << EOF
reset
set title "Eigenstates of the Quartic Oscillator"
set xlabel "x"
unset key
set xrange [-6:6]
set terminal pdf
set output "q5_wavefunctions_figure.pdf"
plot for [ii=0:9] "q5_output.dat" index ii u 1:2 w l
EOF

echo "Plotting the probability distributions..."
gnuplot --persist << EOF
reset
set title "Probability Density of the Quartic Oscillator"
set xlabel "x"
set ylabel "P(x)"
unset key
set xrange [-6:6]
set terminal pdf
set output "q5_densities_figure.pdf"
plot for [ii=0:9] "q5_output.dat" index ii u 1:3 w l
EOF

echo "Cleaning up..."

sed -ie "8s/option = 1/option = 0/" problem4.f90
sed -ie "31s/file_out_index = 12/file_out_index = 11/" problem4.f90
rm *.f90e

echo "Finished..."
echo "====== Problem 5 END ======"

echo "Removing extraneous files FROM EXISTENCE..."
rm p? p1_? *.dat fort.*

echo "Done. Have a nice day!"