program pde
   use iso_c_binding, only: c_double
   use matrix_utils
   implicit none

   ! Parameters and constants
   integer, parameter :: N = 100, iter_max = 100
   real(c_double), parameter :: L = 1.0d0, alpha = 1.0d-4, dt = 0.1d0
   real(c_double) :: dx, lambda, total_heat
   integer :: i, iter, method_choice
   logical :: use_btcs

   ! Arrays for the heat profile and matrices
   real(c_double), dimension(N) :: u
   real(c_double), dimension(N, N) :: A, B

   ! Interface to the external C routine for random numbers
   interface
      function get_random_double() bind(C, name="get_random_double")
         import :: c_double
         real(c_double) :: get_random_double
      end function
   end interface

   ! Calculate spatial step and lambda stability parameter
   dx = L / N
   lambda = alpha * dt / (dx**2)

   ! Ask the user to choose the method:
   print *, "Enter 1 for FTCS or 2 for BTCS:"
   read(*,*) method_choice
   if (method_choice == 2) then
      use_btcs = .true.
      print *, "BTCS method selected."
   else
      use_btcs = .false.
      print *, "FTCS method selected."
   end if

   ! Assemble the FTCS matrix A (with periodic boundary conditions)
   A = 0.0d0
   do i = 1, N
      A(i, modulo(i-2, N)+1) = lambda
      A(i, i)                 = 1.0d0 - 2.0d0 * lambda
      A(i, modulo(i, N)+1)     = lambda
   end do

   if (use_btcs) then
      ! For BTCS, construct B = 2I - A (i.e. B(i,i)=1+2λ, and off-diagonals = -λ)
      B = 0.0d0
      do i = 1, N
         B(i, i) = 2.0d0
      end do
      B = B - A
      ! Invert matrix B using LAPACK routines
      call invert_matrix(B)
   end if

   ! Initialize the heat profile using random values from the external C routine
   do i = 1, N
      u(i) = get_random_double()
   end do

   ! Compute and display the initial total heat
   total_heat = sum(u) * dx
   print *, "Initial total heat =", total_heat

   ! Time evolution loop
   do iter = 1, iter_max
      if (use_btcs) then
         ! BTCS: u^{n+1} = B^{-1} u^{n}
         u = matmul(B, u)
      else
         ! FTCS: u^{n+1} = A u^{n}
         u = matmul(A, u)
      end if
      total_heat = sum(u) * dx
      print *, "Iteration", iter, "Total heat =", total_heat
   end do

end program pde
