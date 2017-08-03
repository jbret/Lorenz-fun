program lorenz
implicit none

integer, parameter :: rprec = kind(1.d0)
integer :: i = 0, n = 1e5 ! iteration counter, total iterations
real(rprec) :: dt = 0.01  ! time step size
real(rprec), dimension(3) :: x, k1=0, k2=0, k3=0, k4=0, g=0, p=0

! initialization
x = (/ 0.0_rprec, 1.0_rprec, 0.0_rprec /) ! initial speeds of 1.0

open(9, file="out_alt.dat")

! low storage version (requires only three levels of storage)
! see CHQZ87 p. 108  (eq. 4.3.19)
do i = 1, n    ! timestepping

  write(9,*)  dt*(i-1), x

  g = x
  p = rhs(x)
  
  x = x + 1/2.0_rprec * dt * p
  g = p
  p = rhs(x)
  
  x = x + 1/2.0_rprec * dt * (p - g)
  g = 1/6.0_rprec * g
  p = rhs(x) - 1/2.0_rprec * p

  x = x + dt * p
  g = g - p
  p = rhs(x) + 2.0_rprec * p
  
  x = x + dt * (g + 1/6.0_rprec * p)

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

!!$subroutine rhs(x_in, k_out)
!!$implicit none
!!$
!!$integer, parameter :: rprec = kind(1.d0)
!!$real(rprec) :: sigma=10.0, rho=28.0, beta=8.0/3
!!$real(rprec), dimension(3) :: x_in, k_out
!!$
!!$k_out(1) = sigma*(x_in(2) - x_in(1))
!!$k_out(2) = x_in(1)*(rho - x_in(3)) - x_in(2)
!!$k_out(3) = x_in(1)*x_in(2) - beta*x_in(3)
!!$
!!$end subroutine

!!$function rhs(x_in) result(k_out)
!!$implicit none
!!$
!!$integer, parameter :: rprec = kind(1.d0)
!!$real(rprec) :: sigma=10.0, rho=28.0, beta=8.0/3
!!$real(rprec), dimension(3) :: x_in, k_out
!!$
!!$k_out(1) = sigma*(x_in(2) - x_in(1))
!!$k_out(2) = x_in(1)*(rho - x_in(3)) - x_in(2)
!!$k_out(3) = x_in(1)*x_in(2) - beta*x_in(3)
!!$
!!$end function rhs
