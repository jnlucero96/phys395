! gaussj.f90  -  starting point for Gauss-Jordan elimination
! compile with: gfortran -O3 -fdefault-real-8 gaussj.f90

program main
implicit none

integer :: i, j ! declare iterator variables
integer, parameter :: n = 10 ! DO NOT TOUCH THIS

! declare matrices and arrays to be used in gaussj
real, dimension(n,n) :: basis_matrix
real, dimension(n,1) :: x, f_of_x ! arrays that hold x values and f(x) values

! declare matrices and arrays to be used in lsolve as benchmark
real, dimension(n,n) :: basis_matrix2
real, dimension(n,1) :: f_of_x2 ! arrays that hold x values and f(x) values

! set up interval, sampled uniformly from [-1.0, 1.0]
call linspace(x, -1.00, 1.00, n)

do i=1,n
    f_of_x(i,1) = func_of_interest(x(i,1))
    f_of_x2(i,1) = func_of_interest(x(i,1))
end do

do j=1,n
    do i=1,n
        basis_matrix(i,j) = ChebyshevT(x(i,1),j)
        basis_matrix2(i,j) = ChebyshevT(x(i,1),j)
    end do
end do

call gaussj(n, basis_matrix, f_of_x)
call lsolve(n, basis_matrix2, f_of_x2)

do i=1,n
    print *, x(i,1), f_of_x(i,1), f_of_x2(i,1)
end do

contains

subroutine swap(a,b)
    ! subroutine to swap two vectors with one another
    real, dimension(:), intent(inout) :: a,b
    real, dimension(size(a)) :: dummy
    dummy=a
    a=b
    b=dummy
end subroutine swap

subroutine lsolve(n, A, B)
    integer n, pivot(n), status; real A(n,n), B(n)

    ! initialize status to all clear
    status = 0

    ! call appropriate version of xGESV
    select case (kind(A))
        case(4); call sgesv(n, 1, A, n, pivot, B, n, status)
        case(8); call dgesv(n, 1, A, n, pivot, B, n, status)
        case default; call abort
    end select

    ! abort at first sign of trouble
    if (status /= 0) stop "singular matrix in lsolve()"
end subroutine

elemental function func_of_interest(x)
    real :: denominator
    real :: func_of_interest
    real, intent(in) :: x
    denominator = (1.0 + 10.0*x*x)
    func_of_interest = 1.0/denominator
end function func_of_interest

elemental function ChebyshevT(x, n)
    real :: ChebyshevT
    real, intent(in) :: x
    integer, intent(in) :: n
    ChebyshevT = cos(n*acos(x))
end function ChebyshevT

subroutine linspace(array, start_point, stop_point, n)
    real, dimension(n, 1), intent(inout) :: array
    real, intent(in) :: start_point, stop_point
    real :: step
    integer, intent(in) :: n
    integer :: i
    ! set up interval, sampled uniformly from [-1.0, 1.0]
    step = (stop_point-start_point) / (n-1.0)
    do i=1,n
        array(i,1) = start_point + (i-1)*step
    end do
end subroutine linspace

! solve A.x = B using Gauss-Jordan elimination
! A gets destroyed, answer is returned in B
subroutine gaussj(n, a, b)

    integer, intent(in) :: n ! number of elements in b
    integer :: i, j ! iterator variables
    real, dimension(n,n), intent(inout) :: a
    real :: scale_factor
    real :: subtract_factor ! initialize a factor
    real, dimension(n,1), intent(inout) :: b

    integer, dimension(n) :: row, col ! pivot indices
    integer,target :: pivot(2) ! store max pivot index
    ! declare pointers that point to the elements of the pivot vector
    integer, pointer :: pivot_1,pivot_2
    ! initialize the pointers
    pivot_1 => pivot(1)
    pivot_2 => pivot(2)

    forall (i=1:n) row(i) = i
    forall (j=1:n) col(j) = j

    ! ========================================================================
    ! ==========
    ! ==========
    ! ========== BEGIN FORWARD ELIMINATION PHASE
    ! ==========
    ! ==========
    ! ========================================================================

    ! iterate through all of the columns save the last
    do j = 1,n-1

        ! find the pivot in the sub-matrix and store them in the pivot array
        pivot_1 = maxloc(abs(A(j:,j)),1) + (j-1);
        pivot_2 = j;

        ! if pivot_1 /= pivot_2 there needs to be a row swap to get the pivot
        ! on the diagonal
        if (pivot_1 /= pivot_2) then
            call swap(A(pivot_1,:),A(pivot_2,:))
            call swap(B(pivot_1,:),B(pivot_2,:))
        end if

    ! zero out the jth column entries that are under the diagonal
    do i = j+1,n
        scale_factor = (A(i,j)/A(j,j))
        A(i,:) = A(i,:) - scale_factor * A(j,:)
        B(i,1) = B(i,1) - scale_factor * B(j,1)
    end do

    end do ! end of iteration through column

    ! check here to make sure that the matrix is in upper triangular form
    ! by checking the lower triangular part of the matrix to make sure all
    ! the entries there are zero
    do i=1,n
        do j=i+1,n
            ! matrix is not in upper triangular form. Bail out.
            if (A(j,i) >= 1.1920929E-07) call abort
        end do
    end do
    ! ========================================================================
    ! ==========
    ! ==========
    ! ========== BEGIN BACK SUBSTITUTION PHASE
    ! ==========
    ! ==========
    ! ========================================================================

    do i = n,1,-1
        subtract_factor = 0.0
        do j = i+1,n
            subtract_factor = subtract_factor + A(i,j)*B(j,1)
        end do
        B(i,1) = (B(i,1) - subtract_factor)/A(i,i)
    end do

end subroutine

end program main
