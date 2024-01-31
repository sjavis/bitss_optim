SUBROUTINE POTENTIAL(x, e, v, vflag, STEST, RMS, PTEST,BOXTEST)
  implicit none
  logical, intent(in) :: vflag
  double precision, intent(in) :: x(2)
  double precision, intent(out) :: e, v(2), RMS
  LOGICAL STEST, PTEST, BOXTEST
  double precision :: e1, e2, e3, v1(2), v2(2), v3(2)

  integer :: i
  double precision :: peak_args(3,5), args(5), dx, dy, etmp

  peak_args(1,:) = [-3., -1.4, 0., 1., 1.]
  peak_args(2,:) = [-2.,  1.4, 0., 1., 1.]
  peak_args(3,:) = [-1., 0.07, 1., 1., 1.]

  e = 0
  v = [0,0]
  do i = 1, 3
    args = peak_args(i,:)

    dx = (x(1) - args(2)) / args(4)
    dy = (x(2) - args(3)) / args(5)
    etmp = args(1) * exp(- dx**2 - dy**2)
    e = e + etmp

    if (.not. vflag) cycle
    v(1) = v(1) - 2 * dx / args(4) * etmp
    v(2) = v(2) - 2 * dy / args(5) * etmp
  end do
END SUBROUTINE POTENTIAL
