program pde
   use iso_c_binding, only: c_double
   use matrix_utils
   implicit none

   ! Parameters and constants
   integer, parameter :: N = 100, iter_max = 100
   real(c_double), parameter :: L = 1.0d0, alpha = 1.0d-4, dt = 0.55d0
   real(c_double) :: dx, lambda, total_heat
   integer :: i, iter, ftcs_unit, btcs_unit
   character(len=10) :: ftcs_file, btcs_file

   ! Arrays for the heat profile and matrices
   real(c_double), dimension(N) :: u_ftcs, u_btcs, u_initial
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
   print *, "Lambda =", lambda

   ! Set up output filenames
   ftcs_file = "ftcs.csv"
   btcs_file = "btcs.csv"

   ! Open output files
   open(newunit=ftcs_unit, file=ftcs_file, status="replace", action="write")
   open(newunit=btcs_unit, file=btcs_file, status="replace", action="write")
   
   ! Write headers to CSV files
   write(ftcs_unit, '(A)') "Iteration,TotalHeat"
   write(btcs_unit, '(A)') "Iteration,TotalHeat"
   
   ! Assemble the FTCS matrix A (with periodic boundary conditions)
   A = 0.0d0
   do i = 1, N
      A(i, modulo(i-2, N)+1) = lambda
      A(i, i)                 = 1.0d0 - 2.0d0 * lambda
      A(i, modulo(i, N)+1)   = lambda
   end do

   ! For BTCS, construct B = 2I - A
   B = 0.0d0
   do i = 1, N
      B(i, i) = 2.0d0
   end do
   B = B - A
   
   ! Invert matrix B using LAPACK routines
   call invert_matrix(B)

   ! Initialize the heat profile using random values
   do i = 1, N
      u_initial(i) = get_random_double()
   end do
   
   ! Copy initial values to both method arrays
   u_ftcs = u_initial
   u_btcs = u_initial

   ! Compute and display the initial total heat
   total_heat = sum(u_initial) * dx
   print *, "Initial total heat =", total_heat
   write(ftcs_unit, '(I0,",",F18.10)') 0, total_heat
   write(btcs_unit, '(I0,",",F18.10)') 0, total_heat

   ! Time evolution loop
   do iter = 1, iter_max
      ! FTCS: u^{n+1} = A u^{n}
      u_ftcs = matmul(A, u_ftcs)
      total_heat = sum(u_ftcs) * dx
      write(ftcs_unit, '(I0,",",F18.10)') iter, total_heat
      
      ! BTCS: u^{n+1} = B^{-1} u^{n}
      u_btcs = matmul(B, u_btcs)
      total_heat = sum(u_btcs) * dx
      write(btcs_unit, '(I0,",",F18.10)') iter, total_heat
      
      if (mod(iter, 10) == 0) then
         print *, "Iteration", iter, "complete"
      end if
   end do

   ! Close the output files
   close(ftcs_unit)
   close(btcs_unit)
   print *, "Results written to files: ", trim(ftcs_file), " and ", trim(btcs_file)

end program pde