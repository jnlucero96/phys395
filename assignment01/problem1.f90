! gaussj.f90  -  starting point for Gauss-Jordan elimination
! compile with: gfortran -O3 -fdefault-real-8 gaussj.f90

program main
implicit none

integer :: n
real, dimension(3,3) :: A
real, dimension(3,1) :: B

! matrix to invert
! A(1,:) = [1.0, 2.0, 3.0]
! A(2,:) = [3.0, 2.0, 1.0]
! A(3,:) = [0.0, 1.0, 0.0]
! B(:,1)= [1.0, 2.0, 3.0]
!n = size(B) ! number of rows

A(1,:) = [2.0, 5.0, 8.0]
A(2,:) = [4.0, 2.0, 2.0]
A(3,:) = [4.0, 1.0, 2.0]
B(:,1)= [3.0, 9.0, 2.0]
n = size(B) ! number of rows

write(*,*)  A(1,:)
write(*,*)  A(2,:)
write(*,*)  A(3,:)
write(*,*) 

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
    integer :: iteration_number
    real, dimension(n,1), intent(inout) :: b

    integer, dimension(n) :: row, col ! pivot indices
    integer,target :: pivot(2) ! store max pivot index
    integer, pointer :: pivot_1,pivot_2
    pivot_1=>pivot(1)
    pivot_2=>pivot(2)
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
    iteration_number = 1
    do j = 1,n-1

        write(*,*) 'iteration', j
                ! find the pivot in the sub-matrix
        !write(*,*) "max loc", kind(maxloc(abs(A(:,j))));

        write(*,*) "max loc", maxloc(abs(A(iteration_number:,j:)),1);


        pivot_1 = maxloc(abs(A(iteration_number:,j)),1) + (iteration_number-1);
        pivot_2=j;
        write(*,*) "Pivot", pivot

                ! (pivot, j) is the location of the maximum of the column j

                ! swap rows/columns
        if (pivot(1)/=pivot(2)) then

        ! ! Actually perform the swap for the matrix A
        dumc = A(pivot(1),:)
        A(pivot(1),:) = A(pivot(2),:)
        A(pivot(2),:) = dumc

        ! ! Perform the swap for the vector B
        v = B(pivot(1),1)
        B(pivot(1),1) = B(pivot(2),1)
        B(pivot(2),1) = v
        end if

!    row([j,pivot(1)]) = row([pivot(1),j])
!    write(*,*) "row",row
!    col([j,pivot(2)]) = col([pivot(2),j])
!    write(*,*)  "col",col



    do i=1,n
        write(*,*) A(i,:)
    end do
    ! get rid of jth column
    do i = j+1,n
        scale_factor = (A(i,j)/A(j,j))
        write(*,*)  "Scale factor =", scale_factor
        A(i,:) = A(i,:) - scale_factor * A(j,:)

        B(i,1) = B(i,1) - scale_factor * B(j,1)
        write(*,*)  "B in loop:", B
        write(*,*)
    end do

    iteration_number = iteration_number + 1

enddo

    write(*,*) "After elimation"

    ! end of column iteration

        write(*,*)  "Before backsub:"
        do i=1,n
            write(*,*)  A(i,:)
        enddo

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

        write(*,*) "Final B vector", B

end subroutine

end program main
