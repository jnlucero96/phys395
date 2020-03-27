program problem1
implicit none

! plotting the potential gives 3 roots. Need to find those roots
! set up upper and lower bound of guesses
real guess(3)

real x
integer i, j

guess = [-4.5, -2.4, 1.0]

do j=1,3 ! iterate through guesses
    ! initial guess
    x = guess(j)

    write (*,*) "Initial guess:", x
    ! Newton's iteration
    do i = 1,16 ! overkill to ensure proper convergence
        x = x - f(x)/df(x)
    end do
    write (*,*) "Root found:", x
end do

contains

! function to find a root of...
pure function f(x); intent(in) x
    real f, x
    f = cos(x) - x/5.0
end function

! derivative of a function to find a root of
pure function df(x); intent(in) x
    real df, x
    df = (-1.0)*(sin(x) + 0.2)
end function

end program problem1
