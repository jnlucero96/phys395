# NOTES

All code is in the code repository. Run the runme file using "./runme.sh" to
run all of the code

## Question1

Use Newton method to find roots of cos(x)-x/5 == 0 to 12 significant digits. 
Output of roots is straight to stdout.
"

## Question2

Use Newton method to find roots of [d/dx P(x)] since we are interested in the
extrema.
Output is straight to stdout.
"plot_potential_q2_figure.pdf" shows the plot of d/dx P(x)

## Question3

Use gl10 integrator to integrate equations of motion.
Output is to "q3_solution.pdf" showing 3 subplots of position, velocity, and
energy deviation.

## Question4

Use gl10 integrator to integrate equations of motion until sign of velocity
flips indicating a turning point.
Output is to "q4_figure.pdf" that shows how long the simulation goes before 
this turning point is reached as a function of initial energy.

## Question5

Use spectral method to solve the Schrodinger equation 

(d^2/dx^2 + 2V(x)) psi = lambda psi

Here lambda = -2E and so implies E = -0.5*lambda. Output of 10 lowest energy
is straight to stdout. 2 plots generated: "p5_wavefunctions_figure.pdf" shows
the wavefunctions "p5_densities_figure.pdf" shows the densities (ie.
wavefunction**2).