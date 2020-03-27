program problem2
implicit none

! plotting the potential gives 3 roots. Need to find those roots
! set up upper and lower bound of guesses
real guess(3)

real x
integer i, j

guess = [-2.9, -0.25, 1.0]

do j=1,3 ! iterate through guesses
    ! initial guess
    x = guess(j)

    write (*,*) "Initial guess:", x
    ! Newton's iteration
    do i = 1,16 ! overkill to ensure proper convergence
        x = x - f(x)/df(x)
    end do
    write (*,*) "Root found:", x
    write (*,*) "Value of function:", true(x)
end do

contains

pure function true(x); intent(in) x
    real true, x
    true = x*(x*x*x+3*x*x-4*x-3) + 4
end function

! function to find a root of...
pure function f(x); intent(in) x
    real f, x
    f = x*(4*x*x+9*x-8) -3
end function

! derivative of a function to find a root of
pure function df(x); intent(in) x
    real df, x
    df = x*(12*x + 18) -8
end function
end program problem2
