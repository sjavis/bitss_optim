MODULE bfgs_logic
  ! Contains the logic used to drive the L-BFGS algorithm.
  ! This should not need to be modified.
  !
  ! The 'heavy lifting' is done the subroutines in lbfgs_maths,
  ! it is those steps that are open to parallelisation, or the calculation
  ! of the gradient, which is handled in the potential file.
  !
  ! These steps have been modified from the pele code to run in pure fortran,
  ! lbfgs_maths has been used directly.
  USE potential
  implicit none

  ! Result Values
  integer :: n_steps, point, CONV_COUNT, n_grad

  ! Internal Values
  double precision, private, allocatable :: s(:,:), y(:,:), p(:), stp(:)
  double precision, private, allocatable :: x0(:), g0(:)
  double precision, private :: e0, e_initial


  CONTAINS


    SUBROUTINE allocate_quench()
      implicit none

      if (allocated(stp))   deallocate(stp);   allocate(stp(N))
      if (allocated(x0))    deallocate(x0);    allocate(x0(N))
      if (allocated(g0))    deallocate(g0);    allocate(g0(N))
      if (allocated(s))     deallocate(s);     allocate(s(n,m))
      if (allocated(y))     deallocate(y);     allocate(y(n,m))
      if (allocated(p))     deallocate(p);     allocate(p(m))
    END SUBROUTINE allocate_quench


    SUBROUTINE quench()
      ! Run the minimiser.
      ! Keep iterating until we reach preset limit or reach
      ! convergence.

      call allocate_quench()
      lbfgs_iter = 1

      call compute_ev_midlayer()
      rms = norm2(g)/sqrt(dble(n))

      e_initial = e

      do while ((lbfgs_iter <= max_iterations) .and. (.not. is_stop_criteron()))
        call one_iteration()
        lbfgs_iter = lbfgs_iter + 1
      end do
    END SUBROUTINE quench


    SUBROUTINE one_iteration()
      ! Take a single step in the L-BFGS algorithm
      character(len=1024) :: filename
      integer :: ispt, iypt

      if (ramp_stat) call ramp(lbfgs_iter)

      ! Assign copies
      e0 = e
      x0 = x
      g0 = g

      ! Get the search direction
      point = mod(lbfgs_iter-1, m) + 1
      call get_step()
      ! Make sure the step is sensible and take it
      call adjust_step_size()

      ! Update the working arrays
      s(:,point) = x - x0
      y(:,point) = g - g0
      p(point) = 1 / dot_product(y(:,point), s(:,point))

      rms = norm2(g)/sqrt(dble(n))
      ! print*, lbfgs_iter, rms

      ! Print log message
      if (debug .and. (mod(lbfgs_iter,lbfgs_log_iter)==0)) then
        write(filename,'(A,A)') trim(output_dir), '/min.log'
        open(1, file=filename, access='append')
        write(1,'(A,ES15.8,A,ES15.8,A,I7,A,ES10.4)') &
          'Energy   ',e,'   RMS   ', rms, '   after step', lbfgs_iter , '  of size  ', norm2(x-x0)
        close(1)
      end if
    END SUBROUTINE one_iteration


    SUBROUTINE get_step()
      implicit none
      integer :: j1, j2, bound
      double precision :: H0, q(n), a(m), b
      
      if (lbfgs_iter == 1) then
        stp = - H0init * g
      
      else
        bound = min(lbfgs_iter-1, m)

        j2 = merge(point-1, m, point>1)
        H0 = 1 / (p(j2) * sum(y(:,j2)**2))

        q = g
        do j1 = point-1, point-bound, -1
          j2 = mod(j1-1+m, m) + 1
          a(j2) = p(j2) * dot_product(s(:,j2), q)
          q = q - a(j2) * y(:,j2)
        end do
        stp = H0 * q

        do j1 = point-bound, point-1
          j2 = mod(j1-1+m, m) + 1
          b = p(j2) * dot_product(y(:,j2), stp)
          stp = stp + s(:,j2) * (a(j2) - b)
        end do
        stp = - stp

        if (dot_product(stp, g) > 0) stp = -stp;
      end if
    END SUBROUTINE get_step


    SUBROUTINE adjust_step_size()
      ! Make sure that the step is sensible, and take it.
      !
      ! This is known as a backtracking line search
      integer :: n_decrease
      double precision :: step_size

      step_size = norm2(stp)
      if (step_size .gt. max_step_size) stp = stp * max_step_size / step_size

      ! decrease step size until it is accepted
      do n_decrease = 1, 10
        x = x0 + stp
        call compute_ev_midlayer()

        ! If the energy rise is too great then reduce step size.
        if (is_accept_step(e, e0)) then
          return
        else
          stp = stp/3d0
        end if
      end do

      !! If the energy is increasing use gradient descent
      !stp = - H0init * g
      !do n_decrease = 1, 3
      !  x = x0 + stp
      !  call compute_ev_midlayer()

      !  ! If the energy rise is too great then reduce step size.
      !  if (is_accept_step(e, e0)) then
      !    return
      !  else
      !    stp = stp/5d0
      !  end if
      !end do
    END SUBROUTINE adjust_step_size


    logical FUNCTION is_stop_criteron()
      double precision :: e0
      logical :: not_initial
      ! Test to see if convergence has been reached.
      ! At the moment this is done by means of a simple RMS check
      !
      ! N.B. stop_signal == .true. implies convergence.

      ! Check that the state has changed from the initialisation
      not_initial = ((e .lt. e_initial-1d-3*abs(e_initial)) .or. allow_initial .or. lbfgs_iter>max_iterations/10)
      if ((rms .lt. convergence_rms) .and. not_initial .and. lbfgs_iter>1) then
        is_stop_criteron = .true.
      else
        is_stop_criteron = .false.
      end if
    END FUNCTION is_stop_criteron


    logical FUNCTION is_accept_step(E_new, E_old)
      ! Do we accept the new step. By default this is dE < dE_max?
      ! We may change this later on.
      double precision, intent(inout) :: E_old, E_new
      double precision :: dE

      ! Relative or absolute energy check
      if (relative_energy_check .eqv. .true.) then
        if (E_old .eq. 0d0) then
          E_old = 1d-40
        end if
        dE = (E_new - E_old)/abs(E_old)
      else
        dE = E_new - E_old
      end if
      is_accept_step = (dE .lt. dE_max)
    END FUNCTION


END MODULE
