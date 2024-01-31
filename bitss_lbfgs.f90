MODULE bitss_lbfgs

  USE COMMONS, ONLY : nopt
  USE KEY, ONLY : bitsslbfgs_m,  BITSSLBFGS_maxstep
  USE BITSS_POTENTIAL, ONLY : bitss_e, bitss_eg

  IMPLICIT NONE

  INTEGER, PRIVATE :: point, lbfgs_iter
  DOUBLE PRECISION, PRIVATE, ALLOCATABLE :: s(:,:), y(:,:), p(:), stp(:)
  DOUBLE PRECISION, PRIVATE, ALLOCATABLE :: x0(:), g0(:), g(:), x(:)
  DOUBLE PRECISION, PRIVATE :: e, e0, e_initial, rms


  CONTAINS


    SUBROUTINE allocate_quench()
      IF (ALLOCATED(stp))   DEALLOCATE(stp);   ALLOCATE(stp(2*nopt))
      IF (ALLOCATED(x0))    DEALLOCATE(x0);    ALLOCATE(x0(2*nopt))
      IF (ALLOCATED(x))     DEALLOCATE(x);     ALLOCATE(x(2*nopt))
      IF (ALLOCATED(g0))    DEALLOCATE(g0);    ALLOCATE(g0(2*nopt))
      IF (ALLOCATED(g))     DEALLOCATE(g);     ALLOCATE(g(2*nopt))
      IF (ALLOCATED(s))     DEALLOCATE(s);     ALLOCATE(s(2*nopt, 2*BITSSLBFGS_m))
      IF (ALLOCATED(y))     DEALLOCATE(y);     ALLOCATE(y(2*nopt, 2*BITSSLBFGS_m))
      IF (ALLOCATED(p))     DEALLOCATE(p);     ALLOCATE(p(2*BITSSLBFGS_m)) ! rho
    END SUBROUTINE allocate_quench


    FUNCTION minimise(coords)
      USE KEY, ONLY : BITSS_COEFITER
      DOUBLE PRECISION :: COORDS(2*nopt)
      LOGICAL :: minimise
      CALL allocate_quench()

      x = coords
      CALL bitss_eg(x, e, g)
      rms = NORM2(g)/SQRT(DBLE(2*nopt))
      e_initial = e

      lbfgs_iter = 1
      open(10, file="log.txt", position="append")
      DO WHILE ((lbfgs_iter <= BITSS_COEFITER) .and. (.not. check_convergence()))
        ! Assign copies
        e0 = e
        x0 = x
        g0 = g

        point = MOD(lbfgs_iter-1, 2*BITSSLBFGS_m) + 1
        CALL get_step() ! Get the search direction
        CALL adjust_step_size() ! Perform a simple linesearch
        CALL bitss_eg(x, e, g)
        write(10,*) x

        ! Update the working arrays
        s(:,point) = x - x0
        y(:,point) = g - g0
        p(point) = 1 / DOT_PRODUCT(y(:,point), s(:,point))

        rms = NORM2(g)/SQRT(DBLE(2*nopt))
        lbfgs_iter = lbfgs_iter + 1
      END DO
      close(10)

      coords = x
      minimise = check_convergence()
    END FUNCTION minimise


    SUBROUTINE get_step()
      INTEGER :: j1, j2, bound
      DOUBLE PRECISION :: H0, q(2*nopt), a(2*BITSSLBFGS_m), b

      IF (lbfgs_iter == 1) THEN
        H0 = 1d-10
        stp = - H0 * g

      ELSE
        bound = MIN(lbfgs_iter-1, 2*BITSSLBFGS_m)

        j2 = MERGE(point-1, 2*BITSSLBFGS_m, point>1)
        H0 = 1 / (p(j2) * SUM(y(:,j2)**2))

        q = g
        DO j1 = point-1, point-bound, -1
          j2 = MOD(j1-1+2*BITSSLBFGS_m, 2*BITSSLBFGS_m) + 1
          a(j2) = p(j2) * DOT_PRODUCT(s(:,j2), q)
          q = q - a(j2) * y(:,j2)
        END DO
        stp = H0 * q

        DO j1 = point-bound, point-1
          j2 = MOD(j1-1+2*BITSSLBFGS_m, 2*BITSSLBFGS_m) + 1
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
      IF ((BITSSLBFGS_maxstep > 0) .and. (step_size > BITSSLBFGS_maxstep)) THEN
        stp = stp * BITSSLBFGS_maxstep / step_size
      END IF

      ! decrease step size until it is accepted
      DO n_decrease = 1, 10
        x = x0 + stp
        CALL bitss_e(x, e)

        ! If the energy rise is too great then reduce step size.
        IF (accept_step(e, e0)) EXIT
        x = x0 - stp
        CALL bitss_e(x, e)
        stp = stp / 1d0
      END DO
    END SUBROUTINE adjust_step_size


    LOGICAL FUNCTION check_convergence()
      USE KEY, ONLY : BITSSLBFGS_conv
      check_convergence = (rms .lt. BITSSLBFGS_conv)
    END FUNCTION check_convergence


    LOGICAL FUNCTION accept_step(E_new, E_old)
      ! Do we accept the new step? This is dE < dE_max
      USE KEY, ONLY : BITSSLBFGS_DEMAX
      DOUBLE PRECISION, INTENT(IN) :: E_old, E_new
      DOUBLE PRECISION :: dE
      ! LOGICAL :: relative_energy_check

      ! Relative or absolute energy check
      ! IF (relative_energy_check) THEN
      !   IF (E_old <= 0d0) print *, "Warning: Attempting relative energy check with non-positive energy"
      !   dE = (E_new - E_old) / E_old
      ! ELSE
      dE = E_new - E_old
      ! END IF
      accept_step = (dE .lt. BITSSLBFGS_DEMAX)
    END FUNCTION


END MODULE
