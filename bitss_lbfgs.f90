MODULE bitss_lbfgs

  USE COMMONS, ONLY : nopt
  USE KEY, ONLY : bitsslbfgs_m
  USE BITSSMODULE, ONLY : bitss_e, bitss_eg

  IMPLICIT NONE

  INTEGER, PRIVATE :: point
  DOUBLE PRECISION, PRIVATE, ALLOCATABLE :: s(:,:), y(:,:), p(:), stp(:)
  DOUBLE PRECISION, PRIVATE, ALLOCATABLE :: x0(:), g0(:), g(:)
  DOUBLE PRECISION, PRIVATE :: e, e0, e_initial, rms
  ASSOCIATE(n=>2*nopt, m=>bitsslbfgs_m)


  CONTAINS


    SUBROUTINE allocate_quench()
      IF (ALLOCATED(stp))   DEALLOCATE(stp);   ALLOCATE(stp(n))
      IF (ALLOCATED(x0))    DEALLOCATE(x0);    ALLOCATE(x0(n))
      IF (ALLOCATED(g0))    DEALLOCATE(g0);    ALLOCATE(g0(n))
      IF (ALLOCATED(g))     DEALLOCATE(g);     ALLOCATE(g(n))
      IF (ALLOCATED(s))     DEALLOCATE(s);     ALLOCATE(s(n, m))
      IF (ALLOCATED(y))     DEALLOCATE(y);     ALLOCATE(y(n, m))
      IF (ALLOCATED(p))     DEALLOCATE(p);     ALLOCATE(p(m)) ! rho
    END SUBROUTINE allocate_quench


    FUNCTION minimise(coords)
      USE KEY, ONLY : BITSSLBFGS_MAXITER
      DOUBLE PRECISION, INTENT(INOUT) :: COORDS(n)
      LOGICAL, INTENT(OUT) :: minimise
      DOUBLE PRECISION :: e, g(n)
      CALL allocate_quench()

      x = coords
      CALL bitss_eg(x, e, g)
      rms = NORM2(g)/SQRT(DBLE(n))
      e_initial = e

      lbfgs_iter = 1
      DO WHILE ((lbfgs_iter <= BITSSLBFGS_MAXITER) .and. (.not. check_convergence()))
        ! Assign copies
        e0 = e
        x0 = x
        g0 = g

        point = MOD(lbfgs_iter-1, m) + 1
        CALL get_step() ! Get the search direction
        CALL adjust_step_size() ! Perform a simple linesearch
        CALL bitss_eg(x, e, g)

        ! Update the working arrays
        s(:,point) = x - x0
        y(:,point) = g - g0
        p(point) = 1 / DOT_PRODUCT(y(:,point), s(:,point))

        rms = NORM2(g)/SQRT(DBLE(n))
        lbfgs_iter = lbfgs_iter + 1
      END DO
      minimise = check_convergence()
    END FUNCTION minimise


    SUBROUTINE get_step()
      INTEGER :: j1, j2, bound
      DOUBLE PRECISION :: H0, q(n), a(m), b

      IF (lbfgs_iter == 1) THEN
        H0 = 1 / NORM2(g)
        stp = - H0 * g

      ELSE
        bound = MIN(lbfgs_iter-1, m)

        j2 = MERGE(point-1, m, point>1)
        H0 = 1 / (p(j2) * SUM(y(:,j2)**2))

        q = g
        DO j1 = point-1, point-bound, -1
          j2 = MOD(j1-1+m, m) + 1
          a(j2) = p(j2) * DOT_PRODUCT(s(:,j2), q)
          q = q - a(j2) * y(:,j2)
        END DO
        stp = H0 * q

        DO j1 = point-bound, point-1
          j2 = MOD(j1-1+m, m) + 1
          b = p(j2) * DOT_PRODUCT(y(:,j2), stp)
          stp = stp + s(:,j2) * (a(j2) - b)
        END DO
        stp = - stp

        IF (DOT_PRODUCT(stp, g) > 0) stp = -stp;
      END IF
    END SUBROUTINE get_step


    SUBROUTINE adjust_step_size()
      ! A simple backtracking line search
      INTEGER :: n_decrease
      DOUBLE PRECISION :: step_size

      step_size = NORM2(stp)
      IF (step_size .gt. max_step_size) stp = stp * max_step_size / step_size

      ! decrease step size until it is accepted
      DO n_decrease = 1, 10
        x = x0 + stp
        CALL bitss_e(x, e)

        ! If the energy rise is too great then reduce step size.
        IF (accept_step(e, e0)) EXIT
        stp = stp / 10d0
      END DO
    END SUBROUTINE adjust_step_size


    LOGICAL FUNCTION check_convergence()
      USE KEY, ONLY : BITSSLBFGS_CONVERGENCE
      stop_criterion = (rms .lt. BITSSLBFGS_CONVERGENCE)
    END FUNCTION check_convergence


    LOGICAL FUNCTION accept_step(E_new, E_old)
      ! Do we accept the new step? This is dE < dE_max
      USE KEY, ONLY : BITSSLBFGS_DEMAX
      DOUBLE PRECISION, INTENT(IN) :: E_old, E_new
      DOUBLE PRECISION :: dE

      ! Relative or absolute energy check
      IF (relative_energy_check) THEN
        IF (E_old <= 0d0) print *, "Warning: Attempting relative energy check with non-positive energy"
        dE = (E_new - E_old) / E_old
      ELSE
        dE = E_new - E_old
      END IF
      accept_step = (dE .lt. BITSSLBFGS_DEMAX)
    END FUNCTION


END MODULE
