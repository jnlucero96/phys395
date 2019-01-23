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

write(*,*) A(1,:)
write(*,*) A(2,:)
write(*,*) A(3,:)
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
    real, dimension(n,1), intent(inout) :: b
    
    integer, dimension(n) :: row, col ! pivot indices
    integer :: pivot(2) ! store max pivot index
    
    forall (i=1:n) row(i) = i
    forall (j=1:n) col(j) = j
    write(*,*) "row:",row
    write(*,*) "col:",col
    write(*,*)
    
    ! access data like A(row(i),col(j))

    ! find location of maximal value
    
    do j = 1,n-1
        pivot = maxloc(abs(A(j:,j:)))
        write(*,*) "j:",j,"Pivot:",pivot
        write(*,*)
        ! (pivot, j) is the location of the maximum of the column j

        ! index (pivot, :) with first row 
        ! row([j,pivot(1)]) = row([pivot(1),j])

        ! swap rows/columns
        row([j,pivot(1)]) = row([pivot(1),j])
        col([j,pivot(2)]) = col([pivot(2),j])

        write(*,*) "row:", row
        write(*,*) "col:", col
        write(*,*)
    
        ! get rid of first column
        do i = j+1,n
            scale_factor = (A(row(i),col(j))/A(row(j),col(j)))
            A(row(i),:) = A(row(i),:) - scale_factor * A(row(j),:)
            B(row(i),1) = B(row(i),1) - scale_factor * B(row(j),1)
            write(*,*) "B in loop:", B, A(row(i),col(j))/A(row(j),col(j))
            write(*,*)
        end do

        dumc = A(row(j),:)
        A(row(j),:) = A(j,:)
        A(j,:) = dumc

        v = B(row(j),1)
        B(row(j),1) = B(j,1)
        B(j,1) = v
        write(*,*) A(1,:)
        write(*,*) A(2,:)
        write(*,*) A(3,:)
        write(*,*) 
        write(*,*) B
        write(*,*)
    end do

    
    
end subroutine

end program main