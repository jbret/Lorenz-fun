program lorenz
implicit none

integer, parameter :: rprec = kind(1.d0)
integer :: i = 0, n = 1e5 ! iteration counter, total iterations
real(rprec) :: dt = 0.01  ! time step size
real(rprec), dimension(3) :: x, k1=0, k2=0, k3=0, k4=0

! initialization
x = (/ 0.0_rprec, 1.0_rprec, 0.0_rprec /) ! initial speeds of 1.0

open(9, file="out.dat")

do i = 1, n    ! timestepping

  write(9,*)  dt*(i-1), x

  call rhs(x,                     k1)
  call rhs(x + dt/2.0_rprec*k1,   k2)
  call rhs(x + dt/2.0_rprec*k2,   k3)
  call rhs(x + dt*k3,             k4)

  ! RK4
  x = x + dt/6.0_rprec*(k1 + 2._rprec*k2 + 2._rprec*k3 + k4)

enddo

close(9)
print*, "Simulation is complete."

end program

subroutine rhs(x_in, k_out)
implicit none

integer, parameter :: rprec = kind(1.d0)
real(rprec) :: sigma=10.0, rho=28.0, beta=8.0/3
real(rprec), dimension(3) :: x_in, k_out

k_out(1) = sigma*(x_in(2) - x_in(1))
k_out(2) = x_in(1)*(rho - x_in(3)) - x_in(2)
k_out(3) = x_in(1)*x_in(2) - beta*x_in(3)

end subroutine

