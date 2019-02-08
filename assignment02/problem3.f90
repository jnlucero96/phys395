! general linear least squares fit

program problem3
use utilities ! call the module

implicit none

! number of basis functions
integer, parameter :: n = 3

real, dimension(n,n) :: A(n,n)
real, dimension(n) :: B, U, S
real, dimension(:), allocatable :: x_pts, y_pts
real, dimension(:), allocatable :: x_data, y_data, fit_vals, chi_sq

integer :: i, j, counter
integer :: status = 0

! initialize accumulators
A = 0.0
U = 0.0

! initialize the counter
counter = 1

allocate(x_pts(500000), y_pts(500000)) ! set really big arrays to hold information

! read data from standard input
do while (status == 0)
    read (*,*,iostat=status) x_pts(counter), y_pts(counter)
    if (status < 0) exit
    
    ! evaluate basis functions at point x
    call evalb(x_pts(counter), n, B)

    ! accumulate least square problem matrices
    forall (i=1:n,j=1:n) A(i,j) = A(i,j) + B(i)*B(j); U = U + y_pts(counter)*B
    
    counter = counter + 1
end do

counter = counter - 1 ! because you add one more than you actually need

allocate(x_data(counter), y_data(counter), fit_vals(counter))

x_data = x_pts(:counter)
y_data = y_pts(:counter)

deallocate(x_pts, y_pts) ! deallocate the big arrays

! solve for best fit parameters
call svdsolve_with_vals(n, S, A, U, 1.0e-6)

call dump_to_grid(n, counter, x_data, U, fit_vals)

chi_sq = abs(y_data-fit_vals)
write (*,*) "Chi Square is:", sum(chi_sq)

contains

! basis functions we are fitting
subroutine evalb(x, n, B)
    real, intent(in) :: x
    integer, intent(in) :: n
    real, dimension(n), intent(inout) :: B
    integer :: i
    
    forall (i=1:n) B(i) = cos(2*PI*x*(i-1))
end subroutine evalb

! evaluate chi-squared expression
subroutine dump_to_grid(n, n_chi, x_data, U, fit_vals)
    integer, intent(in) :: n 
    integer, intent(in) :: n_chi
    real, dimension(n_chi), intent(in) :: x_data
    real, dimension(n), intent(in) :: U
    real, dimension(n_chi), intent(inout) :: fit_vals
    real, dimension(n) :: basis
    integer :: i

    do i = 1,n_chi
        basis = 0.0
        call evalb(x_data(i), n, basis)
        fit_vals(i) = sum(U*basis)
    end do
end subroutine dump_to_grid

end program problem3
