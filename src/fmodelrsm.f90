subroutine fmodelrsm(zeta, y, m, r, s, cpar, dpar, loglik, prob)
implicit none
integer, parameter :: dp = kind(1.0d0)
real(dp), intent(in) :: zeta
integer, intent(in) :: y(s,m)
integer, intent(in) :: m, r, s
real(dp), intent(in) :: dpar(m), cpar(r-1)
real(dp), intent(out) :: loglik, prob(m,r)
integer :: i, j, k
real(dp) :: z(s,m)

do j = 1, m
	prob(j,1) = 1
	do k = 2, r
		prob(j,k) = exp((k - 1) * zeta - (k - 1) * dpar(j) - sum(cpar(1:(k-1))))
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
