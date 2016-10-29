subroutine fmodelseq(zeta, y, m, r, s, bpar, loglik, prob)
implicit none
integer, parameter :: dp = kind(1.0d0)
integer, intent(in) :: m, r, s
integer, intent(in) :: y(s,m)
integer :: i, j, k
real(dp), intent(in) :: zeta, bpar(m,r-1)
real(dp), intent(out) :: loglik, prob(m,r)
real(dp) :: z(s,m), sprb(r-1)

do j = 1, m
	do k = 1, r-1
		sprb(k) = 1/(1 + exp(-zeta + bpar(j,k)))
	end do 
	prob(j,1) = 1 - sprb(1)
	if (r .gt. 2) then
		do k = 2, r-1
			prob(j,k) = (1 - sprb(k)) * product(sprb(1:(k-1)))
		end do
	end if
	prob(j,r) = product(sprb)
end do
do i = 1, s
	do j = 1, m
		z(i,j) = prob(j, y(i,j) + 1)
	end do
end do
loglik = log(sum(product(z,2)))

end subroutine