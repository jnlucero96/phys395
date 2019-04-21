program main
implicit none

integer :: i, j
! real, dimension(3,3) :: A
! real, dimension(3,1) :: B

integer, parameter :: n = 100

real, dimension(n,1) :: x ! arrays that hold x values and f(x) values


real, dimension(n,n) :: A
real, dimension(n,1) :: B ! arrays that hold x values and f(x) values

! A(1,:) = [1.0, 2.0, 3.0]
! A(2,:) = [3.0, 2.0, 1.0]
! A(3,:) = [0.0, 1.0, 0.0]
! B(:,1)= [1.0, 2.0, 3.0]

! A(1,:) = [2.0, 5.0, 8.0]
! A(2,:) = [4.0, 2.0, 2.0]
! A(3,:) = [4.0, 1.0, 2.0]
! B(:,1) = [3.0, 9.0, 2.0]

! ! set up interval, sampled uniformly from [-1.0, 1.0]
call linspace(x, -1.0, 1.0, n)


do i=1,n
    B(i,1) = func_of_interest(x(i,1))
end do

do j=1,n
    do i=1,n
        A(i,j) = ChebyshevT(x(i,1),j)
    end do
end do

print *, "Before solving:"
do i=1,n
    print *, x(i,1), B(i,1)
enddo
print *,

call gaussj(A,B)

print *, "After solving:"
do i=1,n
    print *, x(i,1), B(i,1)
enddo

contains

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
    ! Linear equation solution by Gauss-Jordan elimination.
    ! A is an N ×N input coefficient matrix. b is an N × M input matrix
    ! containing M right-hand-side vectors. On output, a is replaced by its
    ! matrix inverse, and b is replaced by the corresponding set of solution
    ! vectors.

    integer :: i,l, x ! declare iterator variables
    integer :: n ! variable to keep track of the size of the problem
    real :: scale_factor ! variable that keeps track of the inv. pivot quantity

    ! =========================================================================
    ! ========= ARRAY DECLARATIONS
    ! =========================================================================

    ! the row in which that pivot element was originally located.
    integer, dimension(size(a,1)) :: index_rows 
    ! the column of the ith  pivot element, is the ith column that is reduced
    integer, dimension(size(a,1)) :: index_cols 

    ! declare an array to be used for keeping track of if the matrix is singular
    ! or not
    integer, dimension(size(a,1)) :: ipiv 
    logical, dimension(size(a,1)) :: lpiv
    real, dimension(size(a,1)) :: dumc

    ! declare an array which keeps track of the location of the pivot
    ! which can be pointed to. Ie. The "target" of the pointer
    integer, target :: irc(2)
    ! initialize pointers to the elements of the array holding the 
    ! location of the pivot
    integer, pointer :: irow,icol
    irow => irc(1)
    icol => irc(2)

    ipiv = 0
    ! ensure that the dimensions of a is square and ensure that the first 
    ! dimension of b is such that it the system is solvable
    n=assert_equality(size(a,1),size(a,2),size(b,1))

    do i=1,n

        ! =====================================================================
        ! ==========
        ! ==========
        ! ========== BEGIN FORWARD ELIMINATION
        ! ==========
        ! ==========
        ! =====================================================================

        lpiv = (ipiv == 0)
         ! find the location of the pivot element
        irc=maxloc(abs(a),outerand(lpiv,lpiv))

        ! if the same column is identified to have a pivot at 
        ! two different times then the matrix is singular so catch it
        ipiv(icol)=ipiv(icol)+1
        if (ipiv(icol) > 1) call nrerror("gaussj: singular matrix (1)")

        ! If index_rows(i)  ̸= index_cols(i) there is an implied column 
        ! interchange. Physically change a and b then if there is

        if (irow /= icol) then
            call swap(a(irow,:),a(icol,:))
            call swap(b(irow,:),b(icol,:))
        end if

        ! record the location in which the pivot element was before the row 
        ! swap occured
        index_rows(i)=irow
        index_cols(i)=icol

        ! if after the row swap there is a diagonal element that is zero
        ! then the matrix is singular so catch it
        if (a(icol,icol) == 0.0) call nrerror("gaussj: singular matrix (2)")

        ! calculate what you need to scale the row by and then set the diagonal
        ! element to 1.0
        scale_factor=1.0/a(icol,icol)
        a(icol,icol)=1.0

        ! scale the rows of the matrix a by the appropriate scaling factor
        ! and do the same for b
        a(icol,:)=a(icol,:)*scale_factor
        b(icol,:)=b(icol,:)*scale_factor

        ! store the pivot column in a dummy variable for use later
        dumc=a(:,icol)

        ! zero out all the elements of the pivot column
        a(:,icol)=0.0

        ! reset the diagonal element to the calculated scale factor
        a(icol,icol)=scale_factor

        a(1:icol-1,:)=a(1:icol-1,:)-outerprod(dumc(1:icol-1),a(icol,:))
        b(1:icol-1,:)=b(1:icol-1,:)-outerprod(dumc(1:icol-1),b(icol,:))
        a(icol+1:,:)=a(icol+1:,:)-outerprod(dumc(icol+1:),a(icol,:))
        b(icol+1:,:)=b(icol+1:,:)-outerprod(dumc(icol+1:),b(icol,:))
    end do

    ! =====================================================================
    ! ==========
    ! ==========
    ! ========== BEGIN BACK SUBSTITUTION
    ! ==========
    ! ==========
    ! =====================================================================

    ! Interchanging pairs of columns in the reverse order that the permutation 
    ! was built up.
    do l=n,1,-1
        call swap(a(:,index_rows(l)),a(:,index_cols(l)))
    end do

end subroutine gaussj

end program main
