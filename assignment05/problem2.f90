program problem2
implicit none

! iterations for bisection cycle
integer, parameter :: iterations = 1000
real, parameter :: tar = 0.0
! direction to shoot. +1.0 for postive, -1.0 for negative directions
real, parameter :: neg = 1.0
! declare spatial step size
real, parameter :: dx = 6.2831853071795864769252867665590057684Q0/2**9
! declare which potential to use. 0 for quadratic, 1 for quartic
integer, parameter :: option = 0
real, parameter :: sep = 1.0

! declare initial guess for the wave function
real, dimension(3) :: init=[0.0, 1.0, 0.0]

! declare iterator variables
integer ii
real q

if (option==0) then
    do ii = 0,9
        write(*,*) bisection(ii*1.0, ii*1.0 + 1.0, sep)
    end do
else
    do while (q < 5.0)
        write(*,*) bisection(q*1.0, q*1.0 + 0.5, sep)
        q = q + 0.5
    end do
end if

contains

function bisection(a0, b0, sep); intent(in) a0, b0, sep
    real a0, b0, sep, bisection
    real a, b, c, fa, fb, fc
    real, parameter :: eps=1e-16 ! demand double precision in solution

    ! initial interval
    a = a0; fa = check_boundary(a, sep)
    b = b0; fb = check_boundary(b, sep)

    if (fa*fb > 0.0) stop "The root is not bracketed, bailing out..."

    ! bisect the interval
    do while (abs(b-a) > eps*1.0)
        c = (a+b)/2.0; fc = check_boundary(c, sep); if (fc == 0.0) exit

        ! Ridder's variation on the basic method
        c = c + (c-a)*sign(1.0,fa-fb)*fc/sqrt(fc*fc-fa*fb); fc = check_boundary(c, sep)

        if (fc == 0.0) exit
        if (fa*fc < 0.0) then; b = c; fb = fc; end if
        if (fc*fb < 0.0) then; a = c; fa = fc; end if
    end do

    bisection = c
end function bisection

! check the values of the wavefunction at the boundaries
function check_boundary(E, sep); intent(in) E, sep
    real sep, E, check_boundary
    real y(3), y_neg(3)
    integer l

    y = init
    y_neg = -1*init

    ! for a given energy...
    do l=1,2**10
        call gl10(y, dx, 1.0, E) ! evolve to positive x
        call gl10(y_neg, dx, -1.0, E) ! evolve to negative x
    end do

    ! check the boundary values of the two
    check_boundary = (y(2) + sep*y_neg(2)) - tar
end function check_boundary

! potential
function V(x, option); intent(in) x, option
    real V, x
    integer option

    if (option==0) then
        V = 0.5*(x*x)
    else
        V = 0.25*(x**4)
    end if

end function V

! evaluate derivatives
subroutine evalf(u, dudt, neg, E)
        real u(3), dudt(3), E, neg

        associate( x => u(1), y => u(2), yprime => u(3) )
        dudt(1)=neg*1 ! dt = dx
        ! derive these from shrodinger equation
        dudt(2)=yprime
        dudt(3)=2.0*(V(x, option)-E)*y
        end associate
end subroutine evalf

! 10th order implicit Gauss-Legendre integrator
subroutine gl10(y, dt, neg, E)
        integer, parameter :: s = 5, n = 3
        real y(n), g(n,s), dt, E, neg; integer i, k

        ! Butcher tableau for 10th order Gauss-Legendre method
        real, parameter :: a(s,s) = reshape((/ &
                  0.5923172126404727187856601017997934066Q-1, -1.9570364359076037492643214050884060018Q-2, &
                  1.1254400818642955552716244215090748773Q-2, -0.5593793660812184876817721964475928216Q-2, &
                  1.5881129678659985393652424705934162371Q-3,  1.2815100567004528349616684832951382219Q-1, &
                  1.1965716762484161701032287870890954823Q-1, -2.4592114619642200389318251686004016630Q-2, &
                  1.0318280670683357408953945056355839486Q-2, -2.7689943987696030442826307588795957613Q-3, &
                  1.1377628800422460252874127381536557686Q-1,  2.6000465168064151859240589518757397939Q-1, &
                  1.4222222222222222222222222222222222222Q-1, -2.0690316430958284571760137769754882933Q-2, &
                  4.6871545238699412283907465445931044619Q-3,  1.2123243692686414680141465111883827708Q-1, &
                  2.2899605457899987661169181236146325697Q-1,  3.0903655906408664483376269613044846112Q-1, &
                  1.1965716762484161701032287870890954823Q-1, -0.9687563141950739739034827969555140871Q-2, &
                  1.1687532956022854521776677788936526508Q-1,  2.4490812891049541889746347938229502468Q-1, &
                  2.7319004362580148889172820022935369566Q-1,  2.5888469960875927151328897146870315648Q-1, &
                  0.5923172126404727187856601017997934066Q-1/), (/s,s/))
        real, parameter ::   b(s) = [ &
                  1.1846344252809454375713202035995868132Q-1,  2.3931433524968323402064575741781909646Q-1, &
                  2.8444444444444444444444444444444444444Q-1,  2.3931433524968323402064575741781909646Q-1, &
                  1.1846344252809454375713202035995868132Q-1]

        ! iterate trial steps
        g = 0.0; do k = 1,16
                g = matmul(g,a)
                do i = 1,s
                        call evalf(y + g(:,i)*dt, g(:,i), neg, E)
                end do
        end do

        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine gl10

! integrate equations of motion for a given amount of time
function integrate(E, t, dt, neg)
    intent(in) t, dt, neg, E
    real t, dt, integrate, neg, E
    real u(3); integer i, n

    ! vector holds [x, psi, psi prime] initial guess
    u = [0.0, 1.0, 0.0]; ! even solutions
    ! u = [0.0, 0.0, 1.0]; ! odd solutions

    ! number of time steps needed
    n = floor(t/dt)

    do i = 1,n
        call gl10(u, dt, neg, E); if (abs(u(2)) > 100.0) exit
        ! if (mod(i, 100) == 0) write(*,*) u(1), u(2), u(3)
    end do

    call gl10(u,t-n*dt, neg, E)

    ! return wavefunction at time t
    integrate = u(2)
end function integrate

end program problem2
