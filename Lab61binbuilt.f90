program FTCS
   implicit none
   integer, parameter :: dp = selected_real_kind(15, 300)
   integer, parameter :: N = 100, iter_max = 100
   real(kind=dp), parameter :: L = 1.0_dp, alpha = 1.0e-4_dp, dt = 0.1_dp
   real(kind=dp) :: dx, lambda, total_heat
   integer :: i, iter
   real(kind=dp), dimension(N) :: u
   real(kind=dp), dimension(N, N) :: A

   dx = L / N
   lambda = alpha * dt / (dx**2)

   A = 0.0_dp
   do i = 1, N
     A(i, modulo(i-2, N)+1) = lambda
     A(i, i)                = 1.0_dp - 2.0_dp * lambda
     A(i, modulo(i, N)+1)    = lambda
   end do

   ! Initialise the heat profile with random values 
   call random_number(u)

   total_heat = sum(u) * dx
   print *, "Initial total heat =", total_heat

   ! Evolve the system
   do iter = 1, iter_max
      u = matmul(A, u)
      total_heat = sum(u) * dx
      print *, "Iteration", iter, "Total heat =", total_heat
   end do

end program FTCS
