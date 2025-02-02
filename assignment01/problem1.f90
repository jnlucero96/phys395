! gaussj.f90  -  starting point for Gauss-Jordan elimination
! compile with: gfortran -O3 -fdefault-real-8 gaussj.f90

program main
! program that solves the matrix equation Ax=b.
! program outputs data in two columns. First column is the solution from the
! algorithm I designed. The second column is the solution that is
! given by the lsolve algorithm from the lapack package. If the algorithm
! that I coded is correct then it should be the case that they are the same.
implicit none

integer :: n, i
real, dimension(3,3) :: A, A2
real, dimension(3,1) :: B, B2

! matrix to invert
A(1,:) = [1.0, 2.0, 3.0]
A(2,:) = [3.0, 2.0, 1.0]
A(3,:) = [0.0, 1.0, 0.0]
B(:,1)= [1.0, 2.0, 3.0]
n = size(B) ! number of rows
A2(1,:) = [1.0, 2.0, 3.0]
A2(2,:) = [3.0, 2.0, 1.0]
A2(3,:) = [0.0, 1.0, 0.0]
B2(:,1)= [1.0, 2.0, 3.0]

call gaussj(n, A, B)
call lsolve(n, A2, B2)

do i=1,n
    print *, B(i,1), B2(i,1)
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
    integer n, pivot(n), status; real A(n,n), B(n,1)

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
