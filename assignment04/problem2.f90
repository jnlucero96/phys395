program problem2
implicit none

! yes, this is what you think it is
real, parameter :: PI = 3.141592653589793238462643383279502884197169399375

! parameters in the potential
real, parameter :: gbl = 1.0
real, parameter :: mll = 1.0

! temporary variables

real dt
real u_final(4)

! timestep resolving fastest timescale in the scan
dt = 1e-4

u_final = integrate(PI/3.0, -PI/3.0, 100.0, dt)

contains

! potential for 2D particle motion
pure function V(x,y); intent(in) x, y
        real V, x, y
        V = -0.5*mll*gbl*(3.0*cos(x)+cos(y))
end function

! total energy of the dynamical system
pure function E(u); intent(in) u
        real u(4), E, T
        real t1_dot, t2_dot

        associate( x => u(1), y => u(2), p1 => u(3), p2 => u(4) )
                t1_dot = (6.0/(mll))*((2.0*p1-3.0*cos(x-y)*p2) / (16.0-9.0*(cos(x-y)**2)))
                t2_dot = (6.0/(mll))*((8.0*p2-3.0*cos(x-y)*p1) / (16.0-9.0*(cos(x-y)**2)))
                T = (1.0/6.0)*mll*(t2_dot**2 + 4.0*(t1_dot**2) + 3.0*t1_dot*t2_dot*cos(x-y))
                E = T + V(x,y)
        end associate
end function

! evaluate derivatives
subroutine evalf(u, dudt)
        real u(4), dudt(4)

        associate( x => u(1), y => u(2), v => u(3), w => u(4) )
                dudt(1) = (6.0/(mll))*((2.0*v-3.0*cos(x-y)*w) / (16.0-9.0*(cos(x-y)**2)))
                dudt(2) = (6.0/(mll))*((8.0*w-3.0*cos(x-y)*v) / (16.0-9.0*(cos(x-y)**2)))
                dudt(3) = -(0.5*(mll))*(dudt(1)*dudt(2)*sin(x-y)+3.0*(gbl)*sin(x))
                dudt(4) = -(0.5*(mll))*(-dudt(1)*dudt(2)*sin(x-y)+(gbl)*sin(y))
        end associate
end subroutine evalf

! 10th order implicit Gauss-Legendre integrator
subroutine gl10(y, dt)
        integer, parameter :: s = 5, n = 4
        real y(n), g(n,s), dt; integer i, k

        ! Butcher tableau for 8th order Gauss-Legendre method
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
        g = 0.0
        do k = 1,16
                g = matmul(g,a)
                do i = 1,s
                        call evalf(y + g(:,i)*dt, g(:,i))
                end do
        end do

        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine gl10

! integrate equations of motion for a given amount of time
function integrate(x, y, t, dt)
        real x, y, t, dt, integrate(4)
        real u(4), E0; integer i, n, counter

        ! start from a given position at rest
        u = [x, y, 0.0, 0.0]

        ! Calculate initial energy
        E0=E(u)

        ! number of time steps needed
        n = floor(t/dt)

        ! initialize counter for num steps to take before flip
        counter=1

        do i = 1,n
                call gl10(u, dt)
                if (mod(i,100)==0) then ! write out results every 100th time step
                        write (*,'(4g24.16)') i*dt, u(1), u(2), (E(u)-E0)/abs(E0)
                end if
        end do

        call gl10(u, t-n*dt)

        ! return state at time t
        integrate = u
end function

end program problem2
