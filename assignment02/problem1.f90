! general linear least squares fit
! compile with gfortran -O3 -fdefault-real-8 leastsq.f90 -llapack

program problem1
use utilities ! call the module

implicit none

! number of basis functions
integer, parameter :: n = 5

real x, y, A(n,n), B(n), U(n)

integer i, j
integer :: status = 0

! initialize accumulators
A = 0.0
U = 0.0

! read data from standard input
do while (status == 0)
    read (*,*,iostat=status) x, y
    if (status < 0) exit
    
    ! evaluate basis functions at point x
    call evalb(x, n, B)
    
    ! accumulate least square problem matrices
    forall (i=1:n,j=1:n) A(i,j) = A(i,j) + B(i)*B(j); U = U + y*B
end do

! solve for best fit parameters
!call lsolve(n, A, U)
call svdsolve(n, A, U, 1.0e-6)  

! output best fit parameters
do i=1,n
    print *, i, U(i)
end do

contains

! basis functions we are fitting
subroutine evalb(x, n, B)
    integer, intent(in) :: n
    real x, B(n)
    integer, dimension(n) :: alpha 

    forall (i=1:n) alpha(i) = i-1

    B = cos(2*PI*x*alpha)

end subroutine

end program problem1