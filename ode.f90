program general_ode_solver 
    implicit none 

    integer, parameter :: n_eq = 3 
    integer, parameter :: steps = 50000
    real(8), parameter :: dt = 0.001 
    real(8) :: y(n_eq) 
    real(8) :: t


    ! Initial Conditions. 
    y(1) = 1.0 !x  
    y(2) = 1.0 !y  
    y(2) = 1.0 !z
    t = 0.0 

    call integrate_ode(y, t, dt, steps, n_eq, rk4, lorenz_system) 

    print *, 'Integration complete!'
    print *, 'Final State:', y 

contains  

    ! INTEGRATION METHODS 
    ! -------------------

    ! RK4 Method. 

    subroutine rk4(y,t,h,n,derivs)
        implicit none 
        integer, intent(in) :: n 
        real(8), intent(inout) :: y(n) 
        real(8), intent(in) :: t, h 

        interface 
            subroutine derivs(t, y, dydt, n) 
                implicit none 
                integer, intent(in) :: n 
                real(8), intent(in) :: t, y(n) 
                real(8), intent(out) :: dydt(n) 

            end subroutine derivs 
        end interface 

        real(8) :: k1(n), k2(n), k3(n), k4(n)  

        ! k1 = h*f(x,y)
        call derivs(t, y, k1, n) 
        k1 = h*k1

        ! k2 = h*f(t + h/2, y + k1/2) 
        call derivs(t + 0.5*h, y + 0.5*k1, k2, n) 
        k2 = h*k2 

        ! k3 = h*f(t + h/2, y + k2/2) 
        call derivs(t + 0.5*h, y + 0.5*k2, k3, n)  
        k3 = h*k3 

        ! k4 = h*f(t + h, y + k3) 
        call derivs(t + h, y + k3, k4, n)
        k4 = h*k4 

        ! Update solution: y_new = y + (k1 + 2*k2 + 3*k3 * k4)/6 
        y = y + (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0

    end subroutine rk4

  ! GENERAL INTEGRATION DRIVER
    ! ==========================
    
    subroutine integrate_ode(y, t, dt, nsteps, neq, method, derivs)
        implicit none
        integer, intent(in) :: nsteps, neq
        real(8), intent(inout) :: y(neq)
        real(8), intent(inout) :: t
        real(8), intent(in) :: dt
        integer :: i
        
        ! Interface for integration method
        interface
            subroutine method(y, t, h, n, derivs)
                implicit none
                integer, intent(in) :: n
                real(8), intent(inout) :: y(n)
                real(8), intent(in) :: t, h
                interface
                    subroutine derivs(t, y, dydt, n)
                        implicit none
                        integer, intent(in) :: n
                        real(8), intent(in) :: t, y(n)
                        real(8), intent(out) :: dydt(n)
                    end subroutine derivs
                end interface
            end subroutine method
        end interface
        
        ! Interface for derivative function
        interface
            subroutine derivs(t, y, dydt, n)
                implicit none
                integer, intent(in) :: n
                real(8), intent(in) :: t, y(n)
                real(8), intent(out) :: dydt(n)
            end subroutine derivs
        end interface
        
        ! Create output file
        call create_output_file()
        
        ! Main integration loop - CLEAN AND NATURAL
        do i = 1, nsteps
            ! Save current state
            call write_state(t, y, neq)
            
            ! Advance one step using chosen method
            call method(y, t, dt, neq, derivs)
            
            ! Update time
            t = t + dt
        end do
        
        ! Write final state
        call write_state(t, y, neq)
        call close_output_file()
        
    end subroutine integrate_ode
    
    ! FILE I/O UTILITIES
    ! ==================
    
    subroutine create_output_file()
        implicit none
        
        ! Open file for writing (unit 10 is a standard choice)
        ! 'replace' means overwrite if exists
        open(unit=10, file='solution.dat', status='replace')
        
        ! Write header comment
        write(10, '(A)') '# Time X Y Z'
        
    end subroutine create_output_file
    
    subroutine write_state(time, state, n)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: time, state(n)
        
        ! Write time and all state variables
        ! Format: F15.6 means 15 characters total, 6 decimal places
        write(10, '(*(F15.6))') time, state
        
    end subroutine write_state
    
    subroutine close_output_file()
        implicit none
        close(10)
        print *, 'Solution written to solution.dat'
    end subroutine close_output_file
    
    ! EXAMPLE: LORENZ SYSTEM
    ! ======================


    subroutine lorenz_system(t, y, dydt, n) 
        implicit none 
        integer, intent(in) :: n 
        real(8), intent(in) :: t, y(n) 
        real(8), intent(out) :: dydt(n) 

        ! Lorenz parameters 

        real(8), parameter :: sigma = 10.0 
        real(8), parameter :: rho = 28.0
        real(8), parameter :: beta = 8.0/3.0

        ! Lorenz Equations 
        ! dx/dt = σ(x-t), dy/dt = x(ρ-z)-y, dz/dt = xy-βz 

        dydt(1) = sigma*(y(2)-y(1)) 
        dydt(2) = y(1)*(rho - y(3)) - y(2) 
        dydt(3) = y(1)*y(2) - beta*y(3) 

    end subroutine lorenz_system 

end program general_ode_solver
























