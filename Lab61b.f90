program FTCS
   use iso_c_binding, only: c_double
   implicit none
   integer, parameter :: N = 100, iter_max = 100
   real(c_double), parameter :: L = 1.0d0, alpha = 1.0d-4, dt = 0.1d0
   real(c_double) :: dx, lambda, total_heat
   integer :: i, iter
   real(c_double), dimension(N) :: u
   real(c_double), dimension(N, N) :: A
 
   interface
     function get_random_double() bind(C, name="get_random_double")
       import :: c_double
       real(c_double) :: get_random_double
     end function
   end interface
 
   dx = L / N
   lambda = alpha * dt / (dx**2)
 
   A = 0.0d0
   do i = 1, N
     A(i, modulo(i-2, N)+1) = lambda
     A(i, i)                 = 1.0d0 - 2.0d0 * lambda
     A(i, modulo(i, N)+1)     = lambda
   end do
 
   ! Initialise the heat profile with random values
   do i = 1, N
     u(i) = get_random_double()
   end do
 
   total_heat = sum(u) * dx
   print *, "Initial total heat =", total_heat
 
   ! Evolve the system
   do iter = 1, iter_max
      u = matmul(A, u)
      total_heat = sum(u) * dx
      print *, "Iteration", iter, "Total heat =", total_heat
   end do
 
end program FTCS
 
