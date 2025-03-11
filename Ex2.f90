program pde
    use iso_c_binding, only: c_double
    use matrix_utils
    implicit none

    integer, parameter :: N = 100, iter_max = 100000
    real(c_double), parameter :: L = 1.0d0, alpha = 1.0d-4
    real(c_double) :: dt, dx, lambda, total_heat
    integer :: i, iter, ftcs_unit, btcs_unit
    character(len=10) :: ftcs_file, btcs_file

    real(c_double), dimension(N) :: u_ftcs, u_btcs, u_initial, temp
    real(c_double), dimension(N, N) :: A, B

    logical :: periodic
    integer :: bc_choice
    real(c_double) :: VL, VR

    interface
       function get_random_double() bind(C, name="get_random_double")
          import :: c_double
          real(c_double) :: get_random_double
       end function
    end interface

    !------------------------------------------------------------------
    ! Select boundary condition type: 1 for Periodic, 2 for Dirichlet
    print *, "Select boundary condition type:"
    print *, "  1: Periodic (ring)"
    print *, "  2: Dirichlet (bar with fixed boundaries)"
    read(*,*) bc_choice

    if (bc_choice == 1) then
       periodic = .true.
    else if (bc_choice == 2) then
       periodic = .false.
       print *, "Enter left boundary value (VL):"
       read(*,*) VL
       print *, "Enter right boundary value (VR):"
       read(*,*) VR
    else
       print *, "Invalid choice, defaulting to Periodic."
       periodic = .true.
    endif
    !------------------------------------------------------------------

    ! Prompt user for dt
    print *, "Input dt:"
    read(*,*) dt

    dx = L/N
    lambda = alpha * dt / (dx**2)
    print *, "Lambda =", lambda

    ! Set up output filenames
    ftcs_file = "ftcs.csv"
    btcs_file = "btcs.csv"

    open(newunit=ftcs_unit, file=ftcs_file, status="replace", action="write")
    open(newunit=btcs_unit, file=btcs_file, status="replace", action="write")

    ! Write CSV headers
    write(ftcs_unit, '(A)') "Iteration,TotalHeat"
    write(btcs_unit, '(A)') "Iteration,TotalHeat"

    !-------------------------------
    ! Construct matrix A for FTCS
    A = 0.0d0
    if (periodic) then
       ! Periodic: use modulo arithmetic to wrap around.
       do i = 1, N
          A(i, modulo(i-2, N)+1) = lambda
          A(i, i)                = 1.0d0 - 2.0d0 * lambda
          A(i, modulo(i, N)+1)   = lambda
       end do
    else
       ! Dirichlet: no wrap-around. Construct interior rows separately.
       ! First row (left boundary)
       A(1,1) = 1.0d0 - 2.0d0 * lambda
       A(1,2) = lambda
       ! Interior rows
       do i = 2, N-1
          A(i, i-1) = lambda
          A(i, i)   = 1.0d0 - 2.0d0 * lambda
          A(i, i+1) = lambda
       end do
       ! Last row (right boundary)
       A(N, N-1) = lambda
       A(N, N)   = 1.0d0 - 2.0d0 * lambda
    endif
    !-------------------------------

    !-------------------------------
    ! Construct BTCS matrix
    B = 0.0d0
    if (periodic) then
       ! For periodic, BTCS uses B = 2I - A.
       do i = 1, N
          B(i, i) = 2.0d0
       end do
       B = B - A
       call invert_matrix(B)
    else
       ! For Dirichlet, BTCS matrix (B̃) is set up without wrap-around.
       ! First row (left boundary)
       B(1,1) = 1.0d0 + 2.0d0 * lambda
       B(1,2) = -lambda
       ! Interior rows
       do i = 2, N-1
          B(i, i-1) = -lambda
          B(i, i)   = 1.0d0 + 2.0d0 * lambda
          B(i, i+1) = -lambda
       end do
       ! Last row (right boundary)
       B(N, N-1) = -lambda
       B(N, N)   = 1.0d0 + 2.0d0 * lambda
       call invert_matrix(B)
    endif
    !-------------------------------

    ! Initialise the heat profile with random values
    do i = 1, N
       u_initial(i) = get_random_double()
    end do
    if (.not. periodic) then
       ! Enforce fixed boundary conditions in the initial state
       u_initial(1) = VL
       u_initial(N) = VR
    endif
    u_ftcs = u_initial
    u_btcs = u_initial

    total_heat = sum(u_initial) * dx
    print *, "Initial total heat =", total_heat
    write(ftcs_unit, '(I0,",",F20.15)') 0, total_heat
    write(btcs_unit, '(I0,",",F20.15)') 0, total_heat

    !-------------------------------
    ! Time evolution loop
    do iter = 1, iter_max

       ! FTCS update
       if (periodic) then
          u_ftcs = matmul(A, u_ftcs)
       else
          u_ftcs = matmul(A, u_ftcs)
          ! Add fixed boundary source terms: λ * VL and λ * VR at the left/right.
          u_ftcs(1) = u_ftcs(1) + lambda * VL
          u_ftcs(N) = u_ftcs(N) + lambda * VR
          ! Enforce the fixed boundary values explicitly.
          u_ftcs(1) = VL
          u_ftcs(N) = VR
       endif
       total_heat = sum(u_ftcs) * dx
       write(ftcs_unit, '(I0,",",F20.15)') iter, total_heat

       ! BTCS update
       if (periodic) then
          u_btcs = matmul(B, u_btcs)
       else
          temp = u_btcs
          ! Add the boundary source terms to the right-hand side vector.
          temp(1) = temp(1) + lambda * VL
          temp(N) = temp(N) + lambda * VR
          u_btcs = matmul(B, temp)
          ! Enforce fixed boundary values.
          u_btcs(1) = VL
          u_btcs(N) = VR
       endif
       total_heat = sum(u_btcs) * dx
       write(btcs_unit, '(I0,",",F20.15)') iter, total_heat

    end do
    !-------------------------------

    close(ftcs_unit)
    close(btcs_unit)
    print *, "Results written to files:", trim(ftcs_file), " and ", trim(btcs_file)

end program pde
