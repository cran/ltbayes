subroutine fmodel4pp(zeta, y, m, s, apar, bpar, cpar, dpar, loglik, prob)
implicit none				
integer, parameter :: dp = kind(1.0d0)
real(dp), parameter :: sqrt2 = sqrt(2.0_dp)	
real(dp), intent(in) :: zeta, apar(m), bpar(m), cpar(m), dpar(m)											
real(dp), intent(out) :: loglik, prob(m,2) 
integer, intent(in) :: m, s, y(s,m)             
integer :: i, j
real(dp) :: z(s,m)

do j = 1, m
 	prob(j,2) = cpar(j) + (dpar(j) - cpar(j)) * (1 + erf(apar(j)*(zeta - bpar(j))/sqrt2))/2
 	prob(j,1) = 1 - prob(j,2)
end do
do i = 1, s
   do j = 1, m
     z(i,j) = prob(j, y(i,j) + 1)
   end do
end do 
loglik = log(sum(product(z, 2)))

end subroutine