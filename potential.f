!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of OPTIM.
!
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!*************************************************************************
!
!  Here we calculate the two-body Morse potential
!
!*************************************************************************
!
! POTENTIAL(COORDS1, E, G1, .TRUE., .FALSE., RMS, .FALSE., .FALSE.
      SUBROUTINE POTENTIAL (COORDS,E,G,GTEST, STEST, RMS, PTEST,BOXTEST)
      IMPLICIT NONE
      USE COMMONS, ONLY: NOPT
      INTEGER  J1, J2, J3, J4, N
      LOGICAL GTEST,STEST, PTEST, BOXTEST
      DOUBLE PRECISION COORDS(NOPT), E, G(NOPT), TEMP, DIST,
     1                 VNEW(NOPT), R, RR(NOPT,NOPT), ERMR(NOPT,NOPT), ERMRM(NOPT,NOPT), ERMRT(NOPT,NOPT),
     2                 X(NOPT), RHO, P2, RMS

      P2=0.0D0
      DO J1=1,N
         RR(J1,J1)=0.0D0
         ERMR(J1,J1)=0.0D0
         ERMRM(J1,J1)=0.0D0
         ERMRT(J1,J1)=0.0D0
         J3=3*(J1-1)
         DO J2=J1+1,N
            J4=3*(J2-1)
            DIST=(X(J4+1)-X(J3+1))**2+(X(J4+2)-X(J3+2))**2+(X(J4+3)-X(J3+3))**2
            DIST=DSQRT(DIST)
            TEMP=DEXP(-RHO*(DIST-1.0D0))
            P2=P2+TEMP*(TEMP-2.0)
            IF (GTEST) THEN
               R=RHO*DIST
               RR(J2,J1)=1.0D0/R
               ERMR(J2,J1)=TEMP
               ERMRM(J2,J1)=TEMP-1.0D0
               ERMRT(J2,J1)=2.0D0*TEMP-1.0D0
               RR(J1,J2)=RR(J2,J1)
               ERMR(J1,J2)=ERMR(J2,J1)
               ERMRM(J1,J2)=ERMRM(J2,J1)
               ERMRT(J1,J2)=ERMRT(J2,J1)
            ENDIF
         ENDDO
      ENDDO

      IF (.NOT.GTEST) RETURN
      CALL MG(N, X, VNEW, RHO, RR, ERMR, ERMRM)


      ! CALL MS(N, X, RHO, RR, ERMR, ERMRM, ERMRT)

      RETURN
      END
!
!*************************************************************************
!
!  Subroutine MG calculates the cartesian gradient analytically for the Morse potential.
!
!*************************************************************************
!
      SUBROUTINE MG(N, X, V, RHO, RR, ERMR, ERMRM)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4
      DOUBLE PRECISION X(3*N), RHO, V(3*N), RR(N,N),
     1                 DUMMY1, ERMR(N,N), ERMRM(N,N)
!
!  First calculate the gradient analytically.
!
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY1=0.0D0
            DO J4=1,N
               DUMMY1=DUMMY1-2.0*RHO*RHO*(X(J3)-X(3*(J4-1)+J2))*
     1               ERMR(J4,J1)*ERMRM(J4,J1)*RR(J4,J1)
            ENDDO
            V(J3)=DUMMY1
         ENDDO
      ENDDO

      RETURN
      END
