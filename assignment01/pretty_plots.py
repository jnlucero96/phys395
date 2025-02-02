#!/usr/bin/python
from sys import argv
from numpy import loadtxt, linspace, empty, cos, arccos, arange, sin, sqrt, pi
from matplotlib.pyplot import subplots

def prep_cheb_grid(x_array, n_array):
    # use broadcasting to vectorize the chebyshev grid calculation
    return cos(n_array[None,:]*arccos(x_array[:,None]))

def prep_D_cheb_grid(x_array, n_array):
    # use broadcasting to vectorize the chebyshev grid calculation
    return n_array[None,:]*sin(n_array[None,:]*arccos(x_array[:,None]))/sqrt(1-x_array[:,None]*x_array[:,None])

def prep_cheb_line(i_array,n):
    return cos((pi/n)*(i_array-0.5))

def problem2():
    # evaluate the true function at a resolution much higher than either of
    # the estimates
    x = linspace(-1.0, 1.0, 1000)
    true_function = 1.0/(1+10*x*x)

    # unload the data from the files and put them into arrays
    x_sample_n10, gaussj_n10, lsolve_n10 = loadtxt("problem2_10n.dat", unpack=True)
    x_sample_n100, gaussj_n100, lsolve_n100 = loadtxt("problem2_100n.dat", unpack=True)

    # Re-calculate the chebyshev grids
    cheb_grid_n10 = prep_cheb_grid(x_sample_n10, arange(1,11))
    cheb_grid_n100 = prep_cheb_grid(x_sample_n100, arange(1,101))

    # calculate the estimates for n = 10
    f_gaussj_n10_estimate = cheb_grid_n10.dot(gaussj_n10)
    f_lsolve_n10_estimate = cheb_grid_n10.dot(lsolve_n10)

    # calculate the estimates for n = 100
    f_gaussj_n100_estimate = cheb_grid_n100.dot(gaussj_n100)
    f_lsolve_n100_estimate = cheb_grid_n100.dot(lsolve_n100)

    # begin plotting
    fig, ax = subplots(2,1, figsize=(10,10), sharex=True)
    ax[0].plot(x, true_function, lw=3.0, color="black", label="True f(x)")
    ax[0].plot(x_sample_n10, f_gaussj_n10_estimate, lw=3.0, ls="dashed", color="red", label="gaussj() estimate")
    ax[0].plot(x_sample_n10, f_lsolve_n10_estimate, lw=3.0, ls="dotted", color="blue", label="lsolve() estimate")
    ax[0].set_ylabel("f(x)", fontsize=20)
    ax[0].tick_params(axis="both", labelsize=20)
    ax[0].legend(loc=0,prop={"size":12})
    ax[0].set_title("Estimate comparison for n = 10", fontsize=24)

    ax[1].plot(x, true_function, lw=3.0, color="black", label="True f(x)")
    ax[1].plot(x_sample_n100, f_gaussj_n100_estimate, lw=3.0, ls="dashed", color="red", label="gaussj() estimate")
    ax[1].plot(x_sample_n100, f_lsolve_n100_estimate, lw=3.0, ls="dotted", color="blue", label="lsolve() estimate")
    ax[1].set_xlabel("x", fontsize=20)
    ax[1].set_ylabel("f(x)", fontsize=20)
    ax[1].tick_params(axis="both", labelsize=20)
    ax[1].legend(loc=0,prop={"size":12})
    ax[1].set_title("Estimate comparison for n = 100", fontsize=24)

    fig.tight_layout()
    fig.savefig("problem2.pdf")

def problem3():
    # evaluate the true derivative function at a resolution much higher than
    # either of the estimates
    x = linspace(-0.99, 0.99, 1000)
    true_function = -(20.0*x)/((1+10*x*x)**2)

    # unload the data from the files and put them into arrays
    x_sample_n10, gaussj_n10, lsolve_n10 = loadtxt("problem3_10n.dat", unpack=True)
    x_sample_n100, gaussj_n100, lsolve_n100 = loadtxt("problem3_100n.dat", unpack=True)

    # Re-calculate the chebyshev grids
    cheb_prime_grid_n10 = prep_D_cheb_grid(x_sample_n10, arange(1,11))
    cheb_prime_grid_n100 = prep_D_cheb_grid(x_sample_n100, arange(1,101))

    # calculate the estimates for n = 10
    f_gaussj_n10_estimate = cheb_prime_grid_n10.dot(gaussj_n10)
    f_lsolve_n10_estimate = cheb_prime_grid_n10.dot(lsolve_n10)

    # calculate the estimates for n = 100
    f_gaussj_n100_estimate = cheb_prime_grid_n100.dot(gaussj_n100)
    f_lsolve_n100_estimate = cheb_prime_grid_n100.dot(lsolve_n100)

    # begin plotting
    fig, ax = subplots(2,1, figsize=(10,10), sharex=True)
    ax[0].plot(x, true_function, lw=3.0, color="black", label="True f'(x)")
    ax[0].plot(
        x_sample_n10, f_gaussj_n10_estimate,
        lw=3.0, ls="dashed", color="red", label="gaussj() estimate"
        )
    ax[0].plot(x_sample_n10, f_lsolve_n10_estimate, lw=3.0, ls="dotted", color="blue", label="lsolve() estimate")
    ax[0].set_ylabel("f(x)", fontsize=20)
    ax[0].tick_params(axis="both", labelsize=20)
    ax[0].legend(loc=0,prop={"size":12})
    ax[0].set_title("n = 10", fontsize=24)

    ax[1].plot(x, true_function, lw=3.0, color="black", label="True f'(x)")
    ax[1].plot(x_sample_n100, f_gaussj_n100_estimate, lw=3.0, ls="dashed", color="red", label="gaussj() estimate")
    ax[1].plot(x_sample_n100, f_lsolve_n100_estimate, lw=3.0, ls="dotted", color="blue", label="lsolve() estimate")
    ax[1].set_xlabel("x", fontsize=20)
    ax[1].set_ylabel("f(x)", fontsize=20)
    ax[1].tick_params(axis="both", labelsize=20)
    ax[1].legend(loc=0,prop={"size":12})
    ax[1].set_title("n = 100", fontsize=24)

    fig.tight_layout()
    fig.savefig("problem3.pdf")

def problem5():
    # evaluate the true derivative function at a resolution much higher than
    # either of the estimates
    x = prep_cheb_line(arange(1,1001), 1000)
    x_D = prep_cheb_line(arange(1,1001), 1000)
    true_function = 1.0/(1+10*x*x)
    true_D_function = -(20.0*x_D)/((1+10*x_D*x_D)**2)

    # unload the data from the files and put them into arrays
    x_sample_n10, gaussj_n10, lsolve_n10, gaussj_d_n10, lsolve_d_n10 = loadtxt("problem5_10n.dat", unpack=True)
    x_sample_n100, gaussj_n100, lsolve_n100, gaussj_d_n100, lsolve_d_n100= loadtxt("problem5_100n.dat", unpack=True)

    # Re-calculate the chebyshev grids
    cheb_grid_n10 = prep_cheb_grid(x_sample_n10, arange(1,11))
    cheb_prime_grid_n10 = prep_D_cheb_grid(x_sample_n10, arange(1,11))
    cheb_grid_n100 = prep_cheb_grid(x_sample_n100, arange(1,101))
    cheb_prime_grid_n100 = prep_D_cheb_grid(x_sample_n100, arange(1,101))

    # calculate the estimates for n = 10
    f_gaussj_n10_estimate = cheb_grid_n10.dot(gaussj_n10)
    f_lsolve_n10_estimate = cheb_grid_n10.dot(lsolve_n10)
    f_gaussj_d_n10_estimate = cheb_prime_grid_n10.dot(gaussj_d_n10)
    f_lsolve_d_n10_estimate = cheb_prime_grid_n10.dot(lsolve_d_n10)

    # calculate the estimates for n = 100
    f_gaussj_n100_estimate = cheb_grid_n100.dot(gaussj_n100)
    f_lsolve_n100_estimate = cheb_grid_n100.dot(lsolve_n100)
    f_gaussj_d_n100_estimate = cheb_prime_grid_n100.dot(gaussj_d_n100)
    f_lsolve_d_n100_estimate = cheb_prime_grid_n100.dot(lsolve_d_n100)

    # begin plotting
    fig, ax = subplots(2, 2, figsize=(10,10), sharex='col', sharey='row')
    ax[0, 0].plot(x, true_function, lw=3.0, color="black", label="True f(x)")
    ax[0, 0].plot(
        x_sample_n10, f_gaussj_n10_estimate,
        lw=3.0, ls="dashed", color="red", label="gaussj() estimate"
        )
    ax[0, 0].plot(x_sample_n10, f_lsolve_n10_estimate, lw=3.0, ls="dotted", color="blue", label="lsolve() estimate")
    ax[0, 0].tick_params(axis="both", labelsize=20)
    # ax[0, 0].legend(loc=0,prop={"size":12})
    ax[0, 0].set_ylabel("f(x)", fontsize=20)
    ax[0, 0].set_title("n = 10", fontsize=24)

    ax[0, 1].plot(x, true_function, lw=3.0, color="black", label="True f(x)")
    ax[0, 1].plot(
        x_sample_n100, f_gaussj_n100_estimate,
        lw=3.0, ls="dashed", color="red", label="gaussj() estimate"
        )
    ax[0, 1].plot(x_sample_n100, f_lsolve_n100_estimate, lw=3.0, ls="dotted", color="blue", label="lsolve() estimate")
    ax[0, 1].tick_params(axis="both", labelsize=20)
    ax[0, 1].legend(loc=1, prop={"size":12})
    ax[0, 1].set_title("n = 100", fontsize=24)

    ax[1, 0].plot(x, true_D_function, lw=3.0, color="black", label="True f'(x)")
    ax[1, 0].plot(x_sample_n10, f_gaussj_d_n10_estimate, lw=3.0, ls="dashed", color="red", label="gaussj() estimate")
    ax[1, 0].plot(x_sample_n10, f_lsolve_d_n10_estimate, lw=3.0, ls="dotted", color="blue", label="lsolve() estimate")
    ax[1, 0].set_ylabel("f'(x)", fontsize=20)
    ax[1, 0].tick_params(axis="both", labelsize=20)
    ax[1, 0].set_xlim([-1.0, 1.0])
    # ax[1, 0].legend(loc=0,prop={"size":12})

    ax[1, 1].plot(x, true_D_function, lw=3.0, color="black", label="True f'(x)")
    ax[1, 1].plot(x_sample_n100, f_gaussj_d_n100_estimate, lw=3.0, ls="dashed", color="red", label="gaussj() estimate")
    ax[1, 1].plot(x_sample_n100, f_lsolve_d_n100_estimate, lw=3.0, ls="dotted", color="blue", label="lsolve() estimate")
    ax[1, 1].set_xlabel("x", fontsize=20)
    ax[1, 1].tick_params(axis="both", labelsize=20)
    ax[1, 1].set_xlim([-1.0, 1.0])
    ax[1, 1].legend(loc=1, prop={"size":12})

    fig.tight_layout()
    fig.savefig("problem5.pdf")
    return 0

def main(problem_to_solve=None):
    # wrapper function to call other functions

    if problem_to_solve.lower() == "problem2":
        problem2()
    elif problem_to_solve.lower() == "problem3":
        problem3()
    elif problem_to_solve.lower() == "problem5":
        problem5()
    else:
        print "Input not understood. Exiting now"
        exit(1)

if __name__ == "__main__":
    problem_to_solve = argv[1]
    main(problem_to_solve)
