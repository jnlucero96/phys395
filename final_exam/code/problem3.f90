program problem3
implicit none

! data array stores the final state vector
real data_array(2)

! temporary variables
real dt

! timestep resolving fastest timescale in the scan
dt = 1e-2

data_array = integrate(20.0, dt)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dynamical system and Gauss-Legendre integrator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! potential for 1D particle motion
pure function V(x); intent(in) x
    real V, x
    
    V = -12.0/(cosh(x)**2)
end function

! total energy of the dynamical system
pure function E(u); intent(in) u
        real u(2), E
        
        associate( alpha => u(1), beta => u(2) )
            E = 0.5*beta*beta + V(alpha)
        end associate
end function

! evaluate derivatives
subroutine evalf(u, dudt)
        real u(2), dudt(2)
        
        associate( alpha => u(1), beta => u(2) )
            dudt(1) = beta
            dudt(2) = (-24.0)*tanh(alpha)/(cosh(alpha)**2)
        end associate
end subroutine evalf

! 10th order implicit Gauss-Legendre integrator
subroutine gl10(y, dt)
        integer, parameter :: s = 5, n = 2
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
        g = 0.0; do k = 1,16
                g = matmul(g,a)
                do i = 1,s
                        call evalf(y + g(:,i)*dt, g(:,i))
                end do
        end do
        
        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine gl10

! integrate equations of motion for a given amount of time
function integrate(t, dt)
    real t, dt, integrate(2)
    real u(2), E0; integer i, n
    
    ! start from a given position at rest
    u = [1.0, 0.0]; E0 = V(1.0)
    
    ! number of time steps needed
    n = floor(t/dt)
    
    do i = 1,n
        call gl10(u, dt)
        ! write out results every 10 steps
        if (mod(i, 10) == 0) write (*,'(4g24.16)') i*dt, u(1), u(2), (E(u)/E0)-1.0
    end do
    
    call gl10(u,t-n*dt)
    
    ! return state at time t
    integrate = u
end function

end program problem3
