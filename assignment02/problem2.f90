! general linear least squares fit
! compile with gfortran -O3 -fdefault-real-8 leastsq.f90 -llapack

program problem2
use utilities ! call the module

implicit none

! number of basis functions
integer, parameter :: n = 3
! number of points on finer grid
integer, parameter :: m = 10000

real :: x, y
real :: x_min, x_max
real, dimension(n,n) :: A(n,n)
real, dimension(n) :: B, U, S

integer :: i, j
integer :: status = 0

! initialize accumulators
A = 0.0
U = 0.0

x_min = 0; x_max = 0

! read data from standard input
do while (status == 0)
    read (*,*,iostat=status) x, y
    if (status < 0) exit

    ! conditions to find the minimum and maximum of the interval 
    if (x < x_min) x_min = x
    if (x > x_max) x_max = x
    
    ! evaluate basis functions at point x
    call evalb(x, n, B)

    ! accumulate least square problem matrices
    forall (i=1:n,j=1:n) A(i,j) = A(i,j) + B(i)*B(j); U = U + y*B
end do

! solve for best fit parameters
CALL svdsolve_with_vals(n, S, A, U, 1.0e-6)

! output best fit parameters
! since S is returned sorted grab the ratio of its largest value to its smallest
print *, "Ratio of largest to smallest singular values:", S(1)/S(n)
print *, "Ratio of largest to smallest best fit parameters:", maxval(U) / minval(U)

open(unit=3, file="polynomial_fit.dat")
    call dump_to_grid(n, x_min, x_max, m, U)
close(3)

contains

! basis functions we are fitting
subroutine evalb(x, n, B)
    integer, intent(in) :: n
    real, intent(in) :: x
    real, dimension(n), intent(inout) :: B
    integer :: i

    forall (i=1:n) B(i) = cos(2*PI*x*(i-1))
end subroutine evalb

! evaluate polynomial expansion
subroutine dump_to_grid(n, x_min, x_max, m, U)
    integer, intent(in) :: n
    integer, intent(in) :: m
    real, intent(in) :: x_min, x_max
    real, dimension(n) :: U
    real, dimension(m) :: x_pts
    real, dimension(n) :: basis
    integer :: i

    call linspace(x_pts, x_min, x_max, m)

    do i = 1,m
        basis = 0.0
        call evalb(x_pts(i), n, basis)
        write (3,*) x_pts(i), sum(U*basis)
    end do
end subroutine dump_to_grid

end program problem2
