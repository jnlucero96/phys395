! general linear least squares fit
! compile with gfortran -O3 -fdefault-real-8 leastsq.f90 -llapack

program problem3
use utilities ! call the module

implicit none

! number of data lines in file
integer, parameter :: file_max = 50000
! number of basis functions
integer, parameter :: n = 10
! number of fine grid points
integer, parameter :: num_grid_fine = 10000

real :: x_min, x_max
real, dimension(n,n) :: A
real, dimension(n) :: B, U
real, dimension(n) :: temp_basis
real, dimension(num_grid_fine) :: fine_grid
real, dimension(:), allocatable :: x_pts, y_pts
real, dimension(:), allocatable :: x_data, y_data, y_data2, fit_vals, fit_vals2

real, dimension(:,:), allocatable :: new_basis_matrix

integer :: i, j, counter, k
integer :: status = 0

! initialize accumulators
A = 0.0
U = 0.0

! initialize the counter
counter = 1

allocate(x_pts(file_max), y_pts(file_max)) ! set really big arrays to hold information

! read data from standard input
do while (status == 0)
    read (*,*,iostat=status) x_pts(counter), y_pts(counter)
    if (status < 0) exit
    
    ! evaluate basis functions at point x
    call evalb(x_pts(counter), n, B)

    ! accumulate least square problem matrices
    forall (i=1:n,j=1:n) A(i,j) = A(i,j) + B(i)*B(j)
    U = U + y_pts(counter)*B
    
    counter = counter + 1
end do

counter = counter - 1 ! because you add one more than you actually need

allocate(x_data(counter), y_data(counter), y_data2(counter))
allocate(fit_vals(num_grid_fine), fit_vals2(num_grid_fine), new_basis_matrix(counter, n))

x_data = x_pts(:counter)
x_min = minval(x_data)
x_max = maxval(x_data)

call linspace(fine_grid, x_min, x_max, num_grid_fine)

y_data = y_pts(:counter)
y_data2 = y_data

deallocate(x_pts, y_pts)

do k=1,counter
    temp_basis = 0.0
    call evalb(x_data(k), n, temp_basis)
    new_basis_matrix(k,:) = temp_basis
end do

! solve for best fit parameters
! method 1
call svdsolve(n, A, U, 1.0e-6)
call dump_expression(n, x_min, x_max, num_grid_fine, U, fit_vals)

open(unit=1, file="finer_grid_svd.dat")
    do i=1,num_grid_fine
        write(1,*) fine_grid(i), fit_vals(i)
    end do
close(1)

! method 2
call lstsq(counter, n, new_basis_matrix, y_data)
call dump_expression(n, x_min, x_max, num_grid_fine, y_data, fit_vals2)

open(unit=2, file="finer_grid_lstsq.dat")
    do i=1,num_grid_fine
        write(2,*) fine_grid(i), fit_vals(i)
    end do
close(2)

contains

! basis functions we are fitting
subroutine evalb(x, n, B)
    real, intent(in) :: x
    integer, intent(in) :: n
    real, dimension(n), intent(inout) :: B
    integer :: i
    
    forall (i=1:n) B(i) = cos(2*PI*x*(i-1))
end subroutine evalb

! minimize |A.x - B|^2 using LAPACK canned SVD routine (A gets destroyed, the answer is returned in B)
! rcond determines the effective rank of A as described in LAPACK docs. Pass -1.0 for machine precision
subroutine lstsq(num_grid_pts, max_order, A, B)
    integer, intent(in) :: num_grid_pts
    integer, intent(in) :: max_order
    real, dimension(num_grid_pts, max_order), intent(inout) :: A 
    real, dimension(num_grid_pts), intent(inout):: B
    integer :: l
    real, dimension(max_order) :: S
    real, dimension(6*num_grid_pts) :: W
    integer :: rank, status
    real, parameter :: rcond = -1.0
    

    l = 3*num_grid_pts + max_order
    
    ! find solution by singular value decomposition
    status = 0; select case (kind(A))
            case(4); call sgelss(num_grid_pts, max_order, 1, A, num_grid_pts, B, num_grid_pts, S, rcond, rank, W, l, status)
            case(8); call dgelss(num_grid_pts, max_order, 1, A, num_grid_pts, B, num_grid_pts, S, rcond, rank, W, l, status)
            case default; call abort
    end select
    
    ! bail at first sign of trouble
    if (status /= 0) call abort
end subroutine lstsq

! evaluate chi-squared expression
subroutine dump_expression(n, x_min, x_max, num_grid_fine, U, fit_vals)
    integer, intent(in) :: n
    integer, intent(in) :: num_grid_fine 
    real, intent(in) :: x_min, x_max
    real, dimension(num_grid_fine), intent(inout) :: fit_vals
    real, dimension(n) :: U
    real, dimension(num_grid_fine) :: x_pts
    real, dimension(n) :: basis
    integer :: i

    call linspace(x_pts, x_min, x_max, num_grid_fine)

    do i = 1,num_grid_fine
        basis = 0.0
        call evalb(x_pts(i), n, basis)
        fit_vals(i) = sum(U*basis)
    end do
end subroutine dump_expression

end program problem3