program lorenz
implicit none

integer, parameter :: rprec = kind(1.d0)
integer :: i = 0, n = 1e5 ! iteration counter, total iterations
real(rprec) :: dt = 0.001  ! time step size
real(rprec), dimension(3) :: x, k1=0, k2=0, k3=0, k4=0, g=0

! initialization
x = (/ 0.0_rprec, 1.0_rprec, 0.0_rprec /) ! initial speeds of 1.0

open(9, file="out_alt3.dat")

! low storage version (requires only two levels of storage)
! see CHQZ87 p. 108-9  (eq. 4.3.20), note typo in last step (forgets the dt...?)
do i = 1, n    ! timestepping

  write(9,*)  dt*(i-1), x

  g = rhs(x)
  
  x = x + 1/3.0_rprec * dt * g
  g = -5/9.0_rprec * g + rhs(x)
  
  x = x + 15/16.0_rprec * dt * g
  g = -153/128.0_rprec * g + rhs(x)

  x = x + 8/15.0_rprec * dt * g

enddo

close(9)
print*, "Simulation is complete."

contains
  function rhs(x_in) result(k_out)
  implicit none

  integer, parameter :: rprec = kind(1.d0)
  real(rprec) :: sigma=10.0, rho=28.0, beta=8.0/3
  real(rprec), dimension(3) :: x_in, k_out

  k_out(1) = sigma*(x_in(2) - x_in(1))
  k_out(2) = x_in(1)*(rho - x_in(3)) - x_in(2)
  k_out(3) = x_in(1)*x_in(2) - beta*x_in(3)

  end function rhs

end program
