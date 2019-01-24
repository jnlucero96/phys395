program main
implicit none

real, dimension(3,3) :: A
real, dimension(3,1) :: B

A(1,:) = [1.0, 2.0, 3.0]
A(2,:) = [3.0, 2.0, 1.0]
A(3,:) = [0.0, 1.0, 0.0]

print *, "Initial Matrix A:"
print *, a(1,:)
print *, a(2,:)
print *, a(3,:)
print *,

B(:,1) = [1.0, 2.0, 3.0]

print *, "Initial Vector B:"
print *, B
print *,

call gaussj(A,B)

print *, "Final Answer:"
print *, B

contains

function assert_equality(n1,n2,n3)
    ! function that checks whether an equality is held
    integer, intent(in) :: n1,n2,n3
    integer :: assert_equality 
    if (n1 == n2 .and. n2 == n3) then
        assert_equality=n1
    else
        stop 'program terminated by assert_equality'
    end if
end function assert_equality 

function outerprod(a,b)
    real, dimension(:), intent(in) :: a,b
    real, dimension(size(a),size(b)) :: outerprod
    outerprod = spread(a,dim=2,ncopies=size(b)) * &
        spread(b,dim=1,ncopies=size(a))
end function outerprod

subroutine nrerror(string)
    character(len=*), intent(in) :: string
    print *, 'nrerror: ',string
    stop 'program terminated by nrerror'
end subroutine nrerror

subroutine swap(a,b)
    ! subroutine to swap two vectors with one another
    real, dimension(:), intent(inout) :: a,b
    real, dimension(size(a)) :: dum
    dum=a
    a=b
    b=dum
end subroutine swap

function outerand(a,b)
    logical, dimension(:), intent(in) :: a,b
    logical, dimension(size(a),size(b)) :: outerand
    outerand = spread(a,dim=2,ncopies=size(b)) .and. &
        spread(b,dim=1,ncopies=size(a))
end function outerand

subroutine gaussj(a,b)

    implicit none

    real, dimension(:,:), intent(inout) :: a, b
    !    Linear equation solution by Gauss-Jordan elimination.
    !    A is an N ×N input coefficient matrix. b is an N × M input matrix
    !    containing M right-hand-side vectors. On output, a is replaced by its
    !    matrix inverse, and b is replaced by the corresponding set of solution
    !    vectors.

    !    These arrays are used for bookkeeping on the pivoting.
    integer, dimension(size(a,1)) :: ipiv,index_rows,index_cols
    logical, dimension(size(a,1)) :: lpiv
    real :: scale_factor
    real, dimension(size(a,1)) :: dumc
    integer, target :: irc(2)
    integer :: i,l,n
    integer, pointer :: irow,icol
    n=assert_equality(size(a,1),size(a,2),size(b,1))

    ! initialize pointers to the elements of the array holding the 
    ! location of the pivot
    irow => irc(1)
    icol => irc(2)

    ipiv=0

    do i=1,n
        lpiv = (ipiv == 0)
         ! find the location of the pivot
        irc=maxloc(abs(a),outerand(lpiv,lpiv))
        ipiv(icol)=ipiv(icol)+1
        if (ipiv(icol) > 1) call nrerror("gaussj: singular matrix (1)")
    !   We now have the pivot element, so we interchange rows, if needed, to put the pivot
    !   element on the diagonal. The columns are not physically interchanged, only relabeled:
    !   index_cols(i), the column of the ith pivot element, is the ith column that is reduced, while
    !   index_rows(i) is the row in which that pivot element was originally located. 


    !   If index_rows(i)  ̸= index_cols(i) there is an implied column interchange. 
    !   With this form of bookkeeping, the solution b’s will end up in the 
    !   correct order, and the inverse matrix will be scrambled by columns.
        if (irow /= icol) then
            call swap(a(irow,:),a(icol,:))
            call swap(b(irow,:),b(icol,:))
        end if
        index_rows(i)=irow
        index_cols(i)=icol
        if (a(icol,icol) == 0.0) &
            call nrerror("gaussj: singular matrix (2)")

        ! calculate what you need to scale the row by
        scale_factor=1.0/a(icol,icol)
        a(icol,icol)=1.0
        a(icol,:)=a(icol,:)*scale_factor
        b(icol,:)=b(icol,:)*scale_factor
        dumc=a(:,icol)
        a(:,icol)=0.0
        a(:,icol)=0.0
        a(icol,icol)=scale_factor
        a(1:icol-1,:)=a(1:icol-1,:)-outerprod(dumc(1:icol-1),a(icol,:))
        b(1:icol-1,:)=b(1:icol-1,:)-outerprod(dumc(1:icol-1),b(icol,:))
        a(icol+1:,:)=a(icol+1:,:)-outerprod(dumc(icol+1:),a(icol,:))
        b(icol+1:,:)=b(icol+1:,:)-outerprod(dumc(icol+1:),b(icol,:))
    end do

    ! It only remains to unscramble the solution in view of the column 
    ! interchanges. We do this by interchanging pairs of columns in the reverse 
    ! order that the permutation was built up.
    do l=n,1,-1
        call swap(a(:,index_rows(l)),a(:,index_cols(l)))
    end do

end subroutine gaussj

end program main
