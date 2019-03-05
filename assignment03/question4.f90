program question4
use utilities
implicit none

contains

pure function df(x, beta, num_pts, num_basis)
    integer, intent(in) :: num_pts, num_basis
    real, intent(in) :: x
    real, dimension(num_pts), intent(in) :: beta
    real, dimension(num_pts,num_basis) :: df
    real :: sum_val
    ! declare iterator variables
    integer :: i, j, ii, jj

    sum_val = 0.0

    do i=1,num_pts
        do j=1,num_basis
            do jj=1,num_basis
                sum_val = sum_val + beta(jj)*basis(x, jj)
            end do
            df(i, j) = basis(x, j)*exp(sum_val)
        end do
    end do
end function df

pure function basis(x, n)
    real, intent(in) :: x
    integer, intent(in) :: n
    real :: basis
    basis = cos(2*PI*n*x)
end function basis

end program question4