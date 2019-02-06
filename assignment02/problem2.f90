! general linear least squares fit
! compile with gfortran -O3 -fdefault-real-8 leastsq.f90 -llapack

program problem2
use utilities ! call the module

implicit none

! number of basis functions
integer, parameter :: n = 100
! number of points on finer grid
integer, parameter :: m = 10000

real :: x, y
real :: x_min, x_max 
real, dimension(9130) :: x_data, y_data
real, dimension(n,n) :: A(n,n) 
real, dimension(n) :: B, U, S

integer :: i, j, counter
integer :: status = 0

! initialize accumulators
A = 0.0
U = 0.0

counter = 1
! read data from standard input
do while (status == 0)
    read (*,*,iostat=status) x, y
    if (status < 0) exit
    
    x_data(counter) = x
    y_data(counter) = y
    
    ! evaluate basis functions at point x
    call evalb(x, n, B)
    
    ! accumulate least square problem matrices
    forall (i=1:n,j=1:n) A(i,j) = A(i,j) + B(i)*B(j); U = U + y*B
end do

x_min = minval(x_data)
x_max = maxval(x_data)

open(unit=1, file="basis_matrix_100n.dat")
    do i=1,n
        write (1,*) A(i,:)
    end do
close(1)

! solve for best fit parameters
!call lsolve(n, A, U)
CALL svdsolve_with_vals(n, S, A, U, 1.0e-6)  

! output best fit parameters
print *, "Calculating the condition number..."
! since S is returned sorted grab the ratio of its largest value to its smallest
print *, "Ratio of largest to smallest singular values:", S(1)/S(n) 
print *, "Ratio of largest to smallest best fit parameters:", maxval(U) / minval(U)

open(unit=2, file="best_fit_params_100n.dat")
    do i=1,n
        write (2,*) U(i)
    end do
close(2)

call dump_expression(n, x_min, x_max, m, U)

contains

! basis functions we are fitting
subroutine evalb(x, n, B)
    integer, intent(in) :: n
    real, intent(in) :: x
    real, dimension(n), intent(inout) :: B
    integer :: i

    forall (i=1:n) B(i) = cos(2*PI*x*i)
end subroutine evalb

! evaluate polynomial expansion
subroutine dump_expression(n, x_min, x_max, m, U)
    integer, intent(in) :: n 
    integer, intent(in) :: m 
    real, intent(in) :: x_min, x_max
    real, dimension(n) :: U
    real, dimension(m,1) :: x_pts
    real, dimension(n) :: basis
    integer :: i 

    call linspace(x_pts, x_min, x_max, m)

    do i = 1,m
        basis = 0.0
        call evalb(x_pts(i,1), n, basis)
        write (*,*) x_pts(i,1), sum(U*basis)
    end do
end subroutine dump_expression

end program problem2