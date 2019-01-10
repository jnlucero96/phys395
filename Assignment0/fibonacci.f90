! author: Joseph Lucero
! date: Thu 12:04:30 10 January 2019
! title: Assignment 0
! compile with: gfortran main.f90 -O3 -o main
! run with: ./main

program fibonacci
implicit none

integer, parameter :: n = 50 ! number of fibonacci numbers to calculate
integer i ! declare iterator variable
integer(8) F(n) ! declare array to store fibonacci numbers in

! Initialize the first two fibonacci numbers
F(1) = 0; F(2) = 1

! write out the first two fibonacci numbers
write(*, *) F(1)
write(*, *) F(2)

! Loop to calculate the rest of the fibonacci numbers
do i = 3,n
    F(i) = F(i-1) + F(i-2) ! fibonacci sequence rule
    write(*, *) F(i) ! write out the fibonacci numbers
enddo

end program fibonacci