subroutine fmodelnrm(zeta, y, m, r, s, apar, bpar, loglik, prob)
implicit none
integer, parameter :: dp = kind(1.0d0)
real(dp), intent(in) :: zeta
integer, intent(in) :: y(s,m), m, r, s
real(dp), intent(in) :: apar(m,r), bpar(m,r)
real(dp), intent(out) :: loglik, prob(m,r)
integer :: i, j, k
real(dp) :: z(s,m)

do j = 1, m
	do k = 1, r
		prob(j,k) = exp(apar(j,k) * zeta + bpar(j,k))
	end do	
	prob(j,:) = prob(j,:)/sum(prob(j,:))
end do
do i = 1, s
	do j = 1, m
		z(i,j) = prob(j, y(i,j) + 1)
	end do
end do
loglik = log(sum(product(z, 2)))

end subroutine
