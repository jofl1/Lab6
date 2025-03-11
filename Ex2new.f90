program pde
    use iso_c_binding, only: c_double
    use matrix_utils
    implicit none

    ! Simulation parameters
    integer, parameter :: N = 100, iter_max = 100
    real(c_double), parameter :: L = 1.0d0, alpha = 1.0d-4
    real(c_double) :: dt, dx, lambda, total_heat
    integer :: ftcs_unit, btcs_unit
    character(len=10) :: ftcs_file = "ftcs.csv", btcs_file = "btcs.csv"

    ! Solution vectors and matrices
    real(c_double), dimension(N) :: u_ftcs, u_btcs, u_initial
    real(c_double), dimension(N, N) :: A, B

    ! Boundary condition variables
    logical :: periodic
    integer :: bc_choice
    real(c_double) :: VL, VR

    ! Interface for external random number generator
    interface
       function get_random_double() bind(C, name="get_random_double")
          import :: c_double
          real(c_double) :: get_random_double
       end function
    end interface

    ! Boundary condition selection
    call select_boundary_conditions(periodic, bc_choice, VL, VR)
    
    ! Temporal discretization
    print *, "Input dt:"
    read(*,*) dt
    
    ! Calculate discretization parameters
    dx = L/N
    lambda = alpha * dt / (dx**2)
    print *, "Lambda =", lambda

    ! Initialize output files
    call initialize_output_files(ftcs_unit, btcs_unit, ftcs_file, btcs_file)

    ! Construct discretization matrices
    call construct_matrices(A, B, lambda, periodic, N)
    
    ! Initialize temperature profile
    call initialize_temperature(u_initial, u_ftcs, u_btcs, periodic, VL, VR, N)
    
    ! Record initial state
    total_heat = sum(u_initial) * dx
    print *, "Initial total heat =", total_heat
    write(ftcs_unit, '(I0,",",F20.15)') 0, total_heat
    write(btcs_unit, '(I0,",",F20.15)') 0, total_heat

    ! Time evolution loop
    call execute_time_evolution(A, B, u_ftcs, u_btcs, dx, lambda, VL, VR, &
                               periodic, N, iter_max, ftcs_unit, btcs_unit)

    ! Cleanup
    close(ftcs_unit)
    close(btcs_unit)
    print *, "Results written to files:", trim(ftcs_file), " and ", trim(btcs_file)

contains

    subroutine select_boundary_conditions(periodic, bc_choice, VL, VR)
        logical, intent(out) :: periodic
        integer, intent(out) :: bc_choice
        real(c_double), intent(out) :: VL, VR
        
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
    end subroutine

    subroutine initialize_output_files(ftcs_unit, btcs_unit, ftcs_file, btcs_file)
        integer, intent(out) :: ftcs_unit, btcs_unit
        character(len=*), intent(in) :: ftcs_file, btcs_file
        
        open(newunit=ftcs_unit, file=ftcs_file, status="replace", action="write")
        open(newunit=btcs_unit, file=btcs_file, status="replace", action="write")
        
        write(ftcs_unit, '(A)') "Iteration,TotalHeat"
        write(btcs_unit, '(A)') "Iteration,TotalHeat"
    end subroutine

    subroutine construct_matrices(A, B, lambda, periodic, N)
        real(c_double), dimension(N,N), intent(out) :: A, B
        real(c_double), intent(in) :: lambda
        logical, intent(in) :: periodic
        integer, intent(in) :: N
        integer :: i
        
        ! Initialize matrices
        A = 0.0d0
        B = 0.0d0
        
        if (periodic) then
            ! FTCS matrix for periodic BC
            do i = 1, N
                A(i, modulo(i-2, N)+1) = lambda
                A(i, i)                = 1.0d0 - 2.0d0 * lambda
                A(i, modulo(i, N)+1)   = lambda
            end do
            
            ! BTCS matrix for periodic BC: B = 2I - A
            do i = 1, N
                B(i, i) = 2.0d0
            end do
            B = B - A
            call invert_matrix(B)
        else
            ! FTCS matrix for Dirichlet BC
            A(1,1) = 1.0d0 - 2.0d0 * lambda
            A(1,2) = lambda
            do i = 2, N-1
                A(i, i-1) = lambda
                A(i, i)   = 1.0d0 - 2.0d0 * lambda
                A(i, i+1) = lambda
            end do
            A(N, N-1) = lambda
            A(N, N)   = 1.0d0 - 2.0d0 * lambda
            
            ! BTCS matrix for Dirichlet BC
            B(1,1) = 1.0d0 + 2.0d0 * lambda
            B(1,2) = -lambda
            do i = 2, N-1
                B(i, i-1) = -lambda
                B(i, i)   = 1.0d0 + 2.0d0 * lambda
                B(i, i+1) = -lambda
            end do
            B(N, N-1) = -lambda
            B(N, N)   = 1.0d0 + 2.0d0 * lambda
            call invert_matrix(B)
        endif
    end subroutine

    subroutine initialize_temperature(u_initial, u_ftcs, u_btcs, periodic, VL, VR, N)
        real(c_double), dimension(N), intent(out) :: u_initial, u_ftcs, u_btcs
        logical, intent(in) :: periodic
        real(c_double), intent(in) :: VL, VR
        integer, intent(in) :: N
        integer :: i
        
        ! Generate random initial temperature profile
        do i = 1, N
            u_initial(i) = get_random_double()
        end do
        
        ! Apply boundary conditions if needed
        if (.not. periodic) then
            u_initial(1) = VL
            u_initial(N) = VR
        endif
        
        ! Initialize solution vectors
        u_ftcs = u_initial
        u_btcs = u_initial
    end subroutine

    subroutine execute_time_evolution(A, B, u_ftcs, u_btcs, dx, lambda, VL, VR, &
                                     periodic, N, iter_max, ftcs_unit, btcs_unit)
        real(c_double), dimension(N,N), intent(in) :: A, B
        real(c_double), dimension(N), intent(inout) :: u_ftcs, u_btcs
        real(c_double), intent(in) :: dx, lambda, VL, VR
        logical, intent(in) :: periodic
        integer, intent(in) :: N, iter_max, ftcs_unit, btcs_unit
        integer :: iter
        real(c_double) :: total_heat
        real(c_double), dimension(N) :: temp
        
        do iter = 1, iter_max
            ! FTCS update
            if (periodic) then
                u_ftcs = matmul(A, u_ftcs)
            else
                u_ftcs = matmul(A, u_ftcs)
                ! Add boundary source terms
                u_ftcs(1) = u_ftcs(1) + lambda * VL
                u_ftcs(N) = u_ftcs(N) + lambda * VR
                ! Enforce boundary values
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
                ! Add boundary source terms
                temp(1) = temp(1) + lambda * VL
                temp(N) = temp(N) + lambda * VR
                u_btcs = matmul(B, temp)
                ! Enforce boundary values
                u_btcs(1) = VL
                u_btcs(N) = VR
            endif
            total_heat = sum(u_btcs) * dx
            write(btcs_unit, '(I0,",",F20.15)') iter, total_heat
        end do
    end subroutine

end program pde