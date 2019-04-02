program problem4
implicit none

! order of the spectral scheme
integer, parameter :: n = 200

! which potential to use. 0 is quadratic, 1 is quartic.
integer, parameter :: option = 0

! scale of compactification
real, parameter :: ell = 1.0

! this is what you think it is...
real, parameter :: pi = 3.1415926535897932384626433832795028842Q0

! collocation grids
real, dimension(n) :: x, theta

! second derivative operator
real, dimension (n,n) :: L, H

real, dimension(n) :: eig_vals_out
real, dimension(n,n) :: eig_vecs_out

integer i, ii, jj
integer file_out_index

eig_vals_out = 0.0
eig_vecs_out = 0.0

file_out_index = 11

! initialize spectral operators
call initg(); call initl()

! Hamiltonian in static Schrödinger equation
H = -L/2.0; forall (i=1:n) H(i,i) = -L(i,i)/2.0 + V(x(i), option)

call lsolve(eig_vals_out, eig_vecs_out)

do ii=1,10
    write(file_out_index,*) eig_vals_out(ii)
    do jj=1,n
        write(*,*) x(jj), eig_vecs_out(jj, ii), eig_vecs_out(jj,ii)**2
    end do
    write(*,*) ""; write(*,*) ""
end do

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! potential in the problem we are trying to solve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! simple harmonic oscillator potential
elemental function V(x, option); intent(in) x, option
    real V, x
    integer option

    if (option==0) then
        V = 0.5*(x*x)
    else
        V = 0.25*(x**4)
    end if

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! spectral grid, basis functions and derivative operators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initialize the collocation grid
subroutine initg()
    integer i

    forall (i=1:n) theta(i) = pi*(n-i+0.5)/n; x = ell/tan(theta)
end subroutine

! evaluate rational Chebyshev basis on a grid theta
subroutine evalb(n, pts, theta, Tn, Tnx, Tnxx)
        integer n, pts; real, dimension(pts), intent(in) :: theta
        real, dimension(pts), intent(out), optional :: Tn, Tnx, Tnxx

        ! Chebyshev basis and its derivatives
        if (present(Tn))   Tn = cos(n*theta)
        if (present(Tnx))  Tnx = n * sin(n*theta) * sin(theta)**2/ell
        if (present(Tnxx)) Tnxx = -n * (n*cos(n*theta)*sin(theta) + 2.0*sin(n*theta)*cos(theta)) * sin(theta)**3/ell**2
end subroutine evalb

! initialize linear spectral derivative operator
subroutine initl()
    integer i, pivot(n), status; real A(n,n), B(n,n)

    ! evaluate basis and differential operator values on collocation grid
    do i = 1,n
        call evalb(i-1, n, theta, Tn=A(i,:), Tnxx=B(i,:))
    end do

        ! find linear operator matrix
        status = 0; select case (kind(A))
                case(4); call sgesv(n, n, A, n, pivot, B, n, status)
                case(8); call dgesv(n, n, A, n, pivot, B, n, status)
                case default; call abort
        end select

        ! bail at first sign of trouble
        if (status /= 0) call abort

        L = transpose(B)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODIFIED Rayleigh itertation solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Rayleigh's iteration solving eigenvalue problem:
! [-1/2 d^2/dx^2 + V(x)] psi(x) = lambda psi(x)
subroutine lsolve(eig_vals, eig_vecs)
    real eig_vals(n), eig_vecs(n,n)
    ! initialize arrays needed for intrinsic lapack function
    real A(n,n), wr(n), wi(n), vl(n,n), vr(n,n), work(4*n)
    integer status

    A = H

    ! find linear operator matrix
    status = 0; select case (kind(A))
            case(4); call sgeev('N','V',n,A,n,wr,wi,vl,1,vr,n,work,4*n,status)
            case(8); call dgeev('N','V',n,A,n,wr,wi,vl,1,vr,n,work,4*n,status)
            case default; call abort
    end select

    ! bail at first sign of trouble
    if (status /= 0) call abort

    ! sort the results according to eigenvalues
    call quicksort(wr, vr, 1, n)

    eig_vals = wr
    eig_vecs = vr
end subroutine lsolve

recursive subroutine quicksort(a, eig_vecs, first, last)
  implicit none
  real  a(n), x, t, eig_vecs(n,n), temp_vec(n)
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     temp_vec = eig_vecs(:,i); eig_vecs(:,i) = eig_vecs(:,j); eig_vecs(:,j) = temp_vec
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, eig_vecs, first, i-1)
  if (j+1 < last)  call quicksort(a, eig_vecs, j+1, last)
end subroutine quicksort

end program problem4
