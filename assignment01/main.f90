! gaussj.f90  -  starting point for Gauss-Jordan elimination
! compile with: gfortran -O3 -fdefault-real-8 gaussj.f90

program main
implicit none

integer :: n
real A(3,3), B(3,1)

! matrix to invert
A(1,:) = [1.0, 2.0, 3.0]
A(2,:) = [3.0, 2.0, 1.0]
A(3,:) = [0.0, 1.0, 0.0]

B(:,1)= [1.0, 2.0, 3.0]
n = size(B) ! number of rows

print *, A(1,:)
print *, A(2,:)
print *, A(3,:)
print *,

call gaussj(n, A, B)

contains

! solve A.x = B using Gauss-Jordan elimination
! A gets destroyed, answer is returned in B
subroutine gaussj(n, a, b)

    integer :: n ! number of elements in b
    integer :: i, j ! iterator variables
    real, dimension(n,n), intent(inout) :: a
    real, dimension(n) :: dumc
    real :: v
    real :: scale_factor
    real :: subtract_factor ! initialize a factor
    real, dimension(n,1), intent(inout) :: b

    integer, dimension(n) :: row, col ! pivot indices
    integer :: pivot(2) ! store max pivot index

    forall (i=1:n) row(i) = i
    forall (j=1:n) col(j) = j
    print *, "row:",row
    print *, "col:",col
    print *,

    ! ========================================================================
    ! ==========
    ! ==========
    ! ========== BEGIN FORWARD ELIMINATION PHASE
    ! ==========
    ! ==========
    ! ========================================================================
    ! iterate through all of the columns save the last
    do j = 1,n-1

        ! find the pivot in the sub-matrix
        pivot = maxloc(abs(A(j:,j:)))
        print *, "j:",j,"Pivot:",pivot
        print *,
        ! (pivot, j) is the location of the maximum of the column j

        ! index (pivot, :) with first row
        ! row([j,pivot(1)]) = row([pivot(1),j])

        ! swap rows/columns
        row([j,pivot(1)]) = row([pivot(1),j])
        col([j,pivot(2)]) = col([pivot(2),j])

        print *, "row:", row
        print *, "col:", col
        print *,

        ! get rid of jth column
        do i = j+1,n
            scale_factor = (A(row(i),col(j))/A(row(j),col(j)))
            print *, "Scale factor =", scale_factor
            A(row(i),:) = A(row(i),:) - scale_factor * A(row(j),:)
            B(row(i),1) = B(row(i),1) - scale_factor * B(row(j),1)
            print *, "B in loop:", B
            print *,
        end do

        ! ! Actually perform the swap for the matrix A
        dumc = A(row(j),:)
        A(row(j),:) = A(j,:)
        A(j,:) = dumc

        ! ! Perform the swap for the vector B
        v = B(row(j),1)
        B(row(j),1) = B(j,1)
        B(j,1) = v

        print *, A(1,:)
        print *, A(2,:)
        print *, A(3,:)
        print *,
        print *, "B outside loop: ", B
        print *,
    end do ! end of column iteration

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

        print *, "Subtract factor:", subtract_factor

        B(i,1) = (B(i,1) - subtract_factor)/A(i,i)
    end do

    print *,"Final B vector", B

end subroutine

end program main