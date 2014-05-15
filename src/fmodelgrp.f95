subroutine fmodelgrp(zeta, y, m, r, s, apar, bpar, loglik, prob)
implicit none
integer, parameter :: dp = kind(1.0d0)
real(dp) :: sqrt2
integer, intent(in) :: m, r, s
real(dp), intent(in) :: zeta, apar(m), bpar(m,r-1)
real(dp), intent(out) :: loglik, prob(m,r)
real(dp) :: z(s,m)
integer, intent(in) :: y(s,m)
integer :: i, j, k

sqrt2 = sqrt(2.0_dp)

do j = 1, m
	prob(j,1) = 1 - (1 + erf(apar(j)*(zeta - bpar(j,1))/sqrt2))/2
	if (r .gt. 2) then
		do k = 2, r-1
			prob(j,k) = 1 - (1 + erf(apar(j)*(zeta - bpar(j,k))/sqrt2))/2 - sum(prob(j,1:(k-1)))
		end do
	end if
	prob(j,r) = 1 - sum(prob(j,1:(r-1)))
end do
do i = 1, s
	do j = 1, m
		z(i,j) = prob(j, y(i,j) + 1)
	end do
end do
loglik = log(sum(product(z,2)))

end subroutine