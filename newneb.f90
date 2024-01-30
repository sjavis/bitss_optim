!   NEB module is an implementation of the nudged elastic band method for performing double-ended pathway searches.
!   Copyright (C) 2003-2006 Semen A. Trygubenko and David J. Wales
!   This file is part of NEB module. NEB module is part of OPTIM.
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
MODULE NEWNEBMODULE
     IMPLICIT NONE
     CONTAINS

     SUBROUTINE NEWNEB(REDOPATH,TSREDO,EINITIAL,QQ,EFINAL,FINFIN,TSRESET)
          USE NEBDATA, ONLY: BADIMAGE, BADPEPTIDE, XYZ, TANVEC, EEE, GGG, SSS, XCART, X, STEPTOT, STARTTIME, SEPARATION, RMS, &
                             NITERDONE, G, ENDTIME, EIMAGE, BADTAU, RRR, XYZCART, GGGCART, XCART,&
                             DVEC, IMGFREEZE, TRUEGRAD, NEWNEBK, DEVIATION, STEPIMAGE
          USE KEYNEB
          USE MINIMISER1
          USE MINIMISER2
! IF_NOT_BIOVIA
          USE PORFUNCS
          USE MINIMISER3
          USE GROWSTRINGUTILS, ONLY: GROWSTRING, TOTSTEPS
          USE GSDATA, ONLY : KEYGSPRINT
          USE MODMEC,ONLY: MECCANOT          
          USE INTCOMMONS, ONLY : DESMINT, INTINTERPT, NINTIM, NDIH, DIHINFO, ALIGNDIR, PREVDIH, NINTC
          USE INTCUTILS, ONLY : INTINTERPOLATE, CART2INT
!          USE AMHGLOBALS, ONLY : NMRES
! hk286
          USE GENRIGID
! END_IF_NOT_BIOVIA
          USE NEBOUTPUT
          USE NEBUTILS
          USE KEY, ONLY : UNRST, GROWSTRINGT, FREEZENODEST, DESMDEBUG, NEBMUPDATE, MUPDATE, BFGSSTEPS, NEBRESEEDT, ORDERI, ORDERJ,&
                          EPSALPHA, REDOBFGSSTEPS, NREPMAX, DISTREF, NEBKINITIAL, ADDREPT, REPPOW, REDOTSIM, MIN1REDO, MIN2REDO,&
                          PUSHOFF, BULKT, D1INIT, D2INIT, REDOPATH1, INTCONSTRAINTT, INTNEBIMAGES, REDOPATH2, PERMDIST, CPPNEBT,&
                          SLERPT, QCIXYZ, QCIDNEBT, NIMAGEMORE, BITSST
!, AMHT, NUMGLY, REDOKADD, RIGIDBODY, TWOD, VARIABLES, WHOLDNEB
          USE MODGUESS,ONLY: GUESSPATHT,NINTERP
          USE SPFUNCTS, ONLY : DUMPCOORDS
          USE NEBTOCONNECT
          USE COMMONS,ONLY: NINTS, PARAM1, PARAM2, PARAM3, REDOPATHNEB, ZSYM, DEBUG, NOPT, NATOMS
          USE MYLBFGS_MOD, ONLY: MYLBFGS
          USE BITSSMODULE ! SJA

          IMPLICIT NONE

          COMMON /OLDC/ EMAX

          DOUBLE PRECISION,INTENT(IN)           :: EINITIAL, EFINAL
          DOUBLE PRECISION                      :: QQ(NOPT),FINFIN(NOPT)

          INTEGER :: J1,JMAX, NPERSIST, ITDONE, K, I, J2, J5, NIMAGENEW, NIMAGEOLD
!, NDONE, NIMAGENEW, NIMAGEOLD
          DOUBLE PRECISION :: EMAX, XDUMMY, TOTALDIST!, LDTOTAL, LDIST, DINCREMENT

          DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: MYPTS ! JMC
          DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: VNEW, LCOORDS ! JMC
          LOGICAL,ALLOCATABLE :: PERSISTENT(:) 
          LOGICAL  PERMDISTSAVE
          LOGICAL REDOPATH, MFLAG, PTEST, LPTEST, LRESET, TSRESET
          DOUBLE PRECISION ENERGY, RMS2, EREAL, TSREDO(*)
!, D, DIST2, RMAT(3,3)
          DOUBLE PRECISION,ALLOCATABLE :: XNEW(:), XOLD(:)

! IF_NOT_BIOVIA
          DOUBLE PRECISION RIGIDQ(DEGFREEDOMS), RIGIDFIN(DEGFREEDOMS) ! sn402
! END_IF_NOT_BIOVIA

          ! efk: for growstrings and internals
          DOUBLE PRECISION, ALLOCATABLE :: DELTAX(:)
          LOGICAL :: GSMFLAG
          LOGICAL :: FAILED

          LOGICAL :: KNOWE, KNOWG, KNOWH 
          COMMON /KNOWN/ KNOWE, KNOWG, KNOWH 
          LOGICAL :: LDEBUG


          IF (BITSST) THEN
            ! SJA> Use BITSS method instead of NEB
            CALL BITSS()
            RETURN
          END IF

          NIMAGEOLD=NIMAGE
          NIMAGENEW=NIMAGE
          IF (READGUESS.AND.(NIMAGE.LT.NIMAGEMORE)) THEN
             NIMAGENEW=NIMAGEMORE
             NIMAGE=NIMAGEMORE ! to get the allocations right
          ENDIF
          
          IF (DESMDEBUG) THEN
          ! output coordinates of endpoints we're trying to connect          
             CALL DUMPCOORDS(QQ,'tryconnect.A.xyz', .FALSE.)
             CALL DUMPCOORDS(FINFIN,'tryconnect.B.xyz', .FALSE.)
          ENDIF
          
          CALL MYCPU_TIME(STARTTIME,.TRUE.)
          IF (NATOMS<=0) THEN
               PRINT '(1x,a)', 'Number of atoms is less or equal to zero. Stop.'
               CALL TSUMMARY
               STOP
          ELSE IF (DEBUG) THEN
               PRINT *, 'newneb> Number of atoms or variables = ',NATOMS
          ENDIF
          ALLOCATE(ORDERI(NREPMAX),ORDERJ(NREPMAX),EPSALPHA(NREPMAX),DISTREF(NREPMAX),REPPOW(NREPMAX))
          ALLOCATE(BADIMAGE(NIMAGE+2),BADPEPTIDE(NIMAGE+2))
          ADDREPT=.FALSE.

          IF (NIMAGE<=0) THEN
               PRINT '(1x,a)', 'Number of images is less or equal to zero. Stop.'
               CALL TSUMMARY
               STOP
          ENDIF
          ! printing
          MOREPRINTING=.FALSE.
          IF (DEBUG.OR.DESMDEBUG) MOREPRINTING=.TRUE.
          IF (MOREPRINTING.AND.(.NOT.INTCONSTRAINTT)) THEN
! IF_NOT_BIOVIA
             IF (GROWSTRINGT) THEN
                CALL KEYGSPRINT(.FALSE.)
! ELSE_IF_BIOVIA
! DOBIOVIA             IF (.FALSE.) THEN
! END_IF_NOT_BIOVIA
             ELSE
                CALL ALLKEYNEBPRINT
             ENDIF
             PRINT*
          ENDIF

          IF (UNRST) THEN
               GRADTYPE="dnebu"
               TANTYPE=4
          ENDIF
          IF (OLDCONNECT) OPTIMIZETS = .FALSE.
          BADTAU=.FALSE.

! set up arrays
          ALLOCATE(XYZ(NOPT*(NIMAGE+2)),GGG(NOPT*(NIMAGE+2)),SSS(NOPT*(NIMAGE+2)),EEE(NIMAGE+2), &
   &               RRR(NIMAGE+2),TANVEC(NOPT,NIMAGE),DVEC(NIMAGE+1),NEWNEBK(NIMAGE+1),DEVIATION(NIMAGE+1),STEPIMAGE(NIMAGE))
          ALLOCATE(PERSISTENT(NIMAGE+2))

            TANVEC(:,:) = 0.0D0
            XYZ(:) = 0.0D0
            EEE(:) = 0.0D0
            GGG(:) = 0.0D0
            SSS(:) = 0.0D0

!         NEWNEBK(1:NIMAGE+1)=NEBK
          NEWNEBK(1:NIMAGE+1)=NEBKINITIAL

! IF_NOT_BIOVIA
          IF (DESMINT) THEN
             ALLOCATE(XYZCART(3*NATOMS*(NIMAGE+2)), GGGCART(3*NATOMS*(NIMAGE+2)), TRUEGRAD(3*NATOMS*(NIMAGE+2)))
             ALLOCATE(DIHINFO(NIMAGE+2,NDIH))
             XCART => XYZCART(3*NATOMS+1:3*NATOMS*(NIMAGE+1))
             DIHINFO(:,:) = 0.0D0
! ELSE_IF_BIOVIA
! DOBIOVIA          IF (.FALSE.) THEN
! END_IF_NOT_BIOVIA
          ELSE
             ALLOCATE(TRUEGRAD(NOPT*(NIMAGE+2)))
          ENDIF          

          X         => XYZ(NOPT+1:NOPT*(NIMAGE+1))
          EIMAGE    => EEE(2:NIMAGE+1)
          G         => GGG(NOPT+1:NOPT*(NIMAGE+1))
          EEE = 0.0D0
          EEE(1)=EINITIAL
          EEE(NIMAGE+2)=EFINAL
! IF_NOT_BIOVIA
          IF (DESMINT.AND..NOT.GROWSTRINGT) THEN
             PRINT*, "newneb>"
             XYZCART(:3*NATOMS) = QQ
             XYZCART(3*NATOMS*(NIMAGE+1)+1:) = FINFIN

             PREVDIH => DIHINFO(1,:)
             CALL CART2INT(QQ,XYZ(:NOPT))

             DO J1 = 2,NIMAGE+2
                ! align all other dihedrals to start
                DIHINFO(J1,:) = DIHINFO(1,:)
             ENDDO

             ALIGNDIR = .TRUE.
             PREVDIH => DIHINFO(NIMAGE+2,:)
             CALL CART2INT(FINFIN,XYZ(NOPT*(NIMAGE+1)+1:))
             ALIGNDIR = .FALSE.
! ELSE_IF_BIOVIA
! DOBIOVIA          IF (.FALSE.) THEN
! END_IF_NOT_BIOVIA
          ELSE
! IF_NOT_BIOVIA
             IF(RIGIDINIT .AND. (RIGIDMOLECULEST .OR. SLERPT)) THEN  ! Should we replace this by IF(RIGIDINIT)?
!
! XYZ contains NIMAGE+2 sets of NOPT coordinates, one for each image
! and each endpoint.
! Fill in the first DEGFREEDOMS coordinates of each image with the
! rigid-body coordinates, then the remainder with 0's (they get ignored
! by calls to the potential anyway)
!
                CALL TRANSFORMCTORIGID(QQ,RIGIDQ)
                CALL TRANSFORMCTORIGID(FINFIN,RIGIDFIN)
                XYZ(:DEGFREEDOMS) = RIGIDQ(:)
                XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+1)+DEGFREEDOMS) = RIGIDFIN(:)
                ATOMRIGIDCOORDT = .FALSE.
! ELSE_IF_BIOVIA
! DOBIOVIA             IF (.FALSE.) THEN
! END_IF_NOT_BIOVIA
             ELSE
                XYZ(1:NOPT)=QQ(1:NOPT)
                XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+2))=FINFIN(1:NOPT)
             ENDIF
          ENDIF
          IF (FREEZENODEST.OR.NEBRESEEDT) THEN
             IF(.NOT.ALLOCATED(IMGFREEZE)) ALLOCATE(IMGFREEZE(NIMAGE))
             IMGFREEZE(:) = .FALSE.
          ENDIF

! IF_NOT_BIOVIA
          IF(GROWSTRINGT) THEN
             IF((DESMINT.AND.NOPT.NE.NINTC).OR.(.NOT.DESMINT.AND.NOPT.NE.3*NATOMS)) THEN
                print*, 'NOPT must be equal to 3*NATOMS or NINTC to use growstring.'
                print*, DESMINT, NINTC, 3*NATOMS, NOPT
                STOP
             ENDIF
             IF (DESMINT) THEN
                CALL GROWSTRING(QQ, FINFIN, NIMAGE, XCART, EIMAGE, TANVEC, RMS, GSMFLAG)
                XYZ(:) = 0.0D0
             ELSE
                CALL GROWSTRING(QQ, FINFIN, NIMAGE, X, EIMAGE, TANVEC, RMS, GSMFLAG)
                CALL IMAGEDISTRIBUTION
             ENDIF
             NITERDONE = TOTSTEPS
! ELSE_IF_BIOVIA
! DOBIOVIA          IF (.FALSE.) THEN
! END_IF_NOT_BIOVIA
          ELSE  ! construct the band

! Read an initial band in from a file, if specified
! IF_NOT_BIOVIA
             IF (READGUESS.OR.(GUESSPATHT.AND.UNRST.AND.(NINTERP.GT.1)).OR.(MECCANOT)) THEN
! ELSE_IF_BIOVIA
! DOBIOVIA             IF (READGUESS.OR.(GUESSPATHT.AND.UNRST.AND.(NINTERP.GT.1))) THEN
! END_IF_NOT_BIOVIA
                NIMAGE=NIMAGEOLD
                CALL RWG("r",.True.,1)
                NIMAGE=NIMAGENEW
                READGUESS=.FALSE.
!
! use the images we have just read in to redistribute NIMAGEMORE images
!
                IF (NIMAGEMORE.GT.NIMAGEOLD) THEN
                   ALLOCATE(XOLD(NOPT*(NIMAGEOLD+2)))
                   XOLD(1:NOPT*(NIMAGEOLD+2))=XYZ(1:NOPT*(NIMAGEOLD+2))
!
! XYZ is reallocated and deallocated in this routine
! DVEC is also deallocated. Do we need to save it? Or recalculate?
!
!                  CALL IMAGEREDISTRIBUTION2(XNEW,NIMAGENEW,XOLD,NIMAGEOLD) 
                   CALL IMAGEREDISTRIBUTION3(XOLD,NIMAGEOLD) 
                   PRINT '(A,I8)','tryconnect> After imageredistribution2 number of interpolation images is ',NIMAGENEW
                   NIMAGE=NIMAGENEW
!                  ALLOCATE(XYZ(NOPT*(NIMAGE+2)),DVEC(NIMAGE+1))
!                  NULLIFY(X)
!                  X         => XYZ(NOPT+1:NOPT*(NIMAGE+1))
!                  XYZ(1:NOPT*(NIMAGE+2))=XNEW(1:NOPT*(NIMAGE+2))
                   DEALLOCATE(XOLD)
                ENDIF

! IF_NOT_BIOVIA
                IF (RIGIDINIT) THEN
                    ! Read coordinates in as atomistic xyz, but we need rigid body coordinates.
                    ! It's possible that the xyz input file doesn't respect the rigid body constraints
                    ! If so, converting to rigid body coordinates will restore the rigid body constraints in the 
                    ! best alignment MINPERMDIST can find to the cartesian input coords.
                    ! This will not work very well if the 'rigid bodies' in the guessed path are very different
                    ! to those in the reference structure.

                    LDEBUG = DEBUG
                    DEBUG = .FALSE. ! We can't have DEBUG set here, because it will detect that the coordinates being
                                    ! read in don't exactly match the rigid body constraints. But this is expected
                                    ! behaviour in this case.
                    CALL GENRIGID_IMAGE_CTORIGID(NIMAGE, XYZ)
                    DEBUG = LDEBUG
                ENDIF 
! END_IF_NOT_BIOVIA
             ELSE
! IF_NOT_BIOVIA
                IF (UNRST) THEN ! JMC
                   ALLOCATE(MYPTS(3*NATOMS*NIMAGE)) ! JMC
                   CALL UNRESDIHENEB(QQ,FINFIN,MYPTS)
                   XYZ(NOPT+1:NOPT*(NIMAGE+1))=MYPTS(1:NOPT*NIMAGE)
                   DEALLOCATE(MYPTS)
! ELSE_IF_BIOVIA
! DOBIOVIA                IF (.FALSE.) THEN
! END_IF_NOT_BIOVIA
                ELSEIF (REDOPATHNEB) THEN
! sn402: There's probably no rigid body support with this option.
!                  REDOKADD=.TRUE.
                   REDOPATH1=.TRUE.
                   REDOTSIM=NIMAGE*D1INIT/(D1INIT+D2INIT)+1
                   XYZ(NOPT*REDOTSIM+1:NOPT*(REDOTSIM+1))=TSREDO(1:NOPT)
                   ALLOCATE(DELTAX(NOPT))
                   IF (BULKT) THEN
                      DO K=1,NATOMS
                         DELTAX(3*(K-1)+1)=XYZ(NOPT*REDOTSIM+3*(K-1)+1) - XYZ(3*(K-1)+1) &
  &                          -PARAM1*NINT((XYZ(NOPT*REDOTSIM+3*(K-1)+1) - XYZ(3*(K-1)+1))/PARAM1)
                         DELTAX(3*(K-1)+2)=XYZ(NOPT*REDOTSIM+3*(K-1)+2) - XYZ(3*(K-1)+2) &
  &                          -PARAM2*NINT((XYZ(NOPT*REDOTSIM+3*(K-1)+2) - XYZ(3*(K-1)+2))/PARAM2)
                         DELTAX(3*(K-1)+3)=XYZ(NOPT*REDOTSIM+3*(K-1)+3) - XYZ(3*(K-1)+3) &
  &                          -PARAM3*NINT((XYZ(NOPT*REDOTSIM+3*(K-1)+3) - XYZ(3*(K-1)+3))/PARAM3)
                      ENDDO
                      DELTAX(1:NOPT)=DELTAX(1:NOPT)/REDOTSIM
                   ELSE
                      DELTAX(1:NOPT) = ( XYZ(NOPT*REDOTSIM+1:NOPT*(REDOTSIM+1) ) - XYZ(1:NOPT) )/REDOTSIM
                   ENDIF
                   DO I=2,REDOTSIM
                      XYZ(NOPT*(I-1)+1:NOPT*I) = XYZ(1:NOPT) + DELTAX*(I-1)
                   ENDDO
                   IF (.NOT.ALLOCATED(VNEW)) ALLOCATE(VNEW(NOPT))
                   IF (.NOT.ALLOCATED(LCOORDS)) ALLOCATE(LCOORDS(NOPT))
                   IF (DEBUG) PRINT '(A)',' newneb> minimising on the start side to make DNEB images'
                   LPTEST=.TRUE.
                   LCOORDS(1:NOPT)=TSREDO(1:NOPT)+PUSHOFF*(MIN1REDO(1:NOPT)-TSREDO(1:NOPT))/D1INIT
                   KNOWE=.FALSE.; KNOWG=.FALSE.
                   LRESET=.TRUE.
                   DO I=REDOTSIM,2,-1
                      CALL MYLBFGS(NOPT,MUPDATE,LCOORDS,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,REDOBFGSSTEPS,LRESET, &
  &                                ITDONE,LPTEST,VNEW,.TRUE.,.FALSE.)
                      LRESET=.FALSE.
                      XYZ(NOPT*(I-1)+1:NOPT*I)=LCOORDS(1:NOPT)
                      IF (DEBUG) PRINT '(A,I6,A,F20.10)',' newneb> energy of image ',I,' is ',EREAL
                   ENDDO
                   REDOPATH1=.FALSE.
                   REDOPATH2=.TRUE.
                   IF (BULKT) THEN
                      DO K=1,NATOMS
                         DELTAX(3*(K-1)+1)=XYZ(NOPT*(NIMAGE+1)+3*(K-1)+1) - XYZ(NOPT*REDOTSIM+3*(K-1)+1) &
  &                          -PARAM1*NINT((XYZ(NOPT*(NIMAGE+1)+3*(K-1)+1) - XYZ(NOPT*REDOTSIM+3*(K-1)+1))/PARAM1)
                         DELTAX(3*(K-1)+2)=XYZ(NOPT*(NIMAGE+1)+3*(K-1)+2) - XYZ(NOPT*REDOTSIM+3*(K-1)+2) &
  &                          -PARAM2*NINT((XYZ(NOPT*(NIMAGE+1)+3*(K-1)+2) - XYZ(NOPT*REDOTSIM+3*(K-1)+2))/PARAM2)
                         DELTAX(3*(K-1)+3)=XYZ(NOPT*(NIMAGE+1)+3*(K-1)+3) - XYZ(NOPT*REDOTSIM+3*(K-1)+3) &
  &                          -PARAM3*NINT((XYZ(NOPT*(NIMAGE+1)+3*(K-1)+3) - XYZ(NOPT*REDOTSIM+3*(K-1)+3))/PARAM3)
                      ENDDO
                      DELTAX(1:NOPT)=DELTAX(1:NOPT)/(NIMAGE-REDOTSIM+1)
                   ELSE
                      DELTAX(1:NOPT) = ( XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+2))    &
  &                                    - XYZ(NOPT*REDOTSIM+1:NOPT*(REDOTSIM+1)) )/(NIMAGE-REDOTSIM+1)
                   ENDIF
                   DO I=REDOTSIM+2,NIMAGE+1
                      XYZ(NOPT*(I-1)+1:NOPT*I) = TSREDO(1:NOPT) + DELTAX*(I-REDOTSIM-1)
                   ENDDO
                   IF (DEBUG) PRINT '(A)',' newneb> minimising on the finish side to make DNEB images'
                   LCOORDS(1:NOPT)=TSREDO(1:NOPT)+PUSHOFF*(MIN2REDO(1:NOPT)-TSREDO(1:NOPT))/D2INIT
                   KNOWE=.FALSE.; KNOWG=.FALSE.
                   LRESET=.TRUE.
                   DO I=REDOTSIM+2,NIMAGE+1
                      CALL MYLBFGS(NOPT,MUPDATE,LCOORDS,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,REDOBFGSSTEPS,LRESET, &
  &                                ITDONE,LPTEST,VNEW,.TRUE.,.FALSE.)
                      LRESET=.FALSE.
                      XYZ(NOPT*(I-1)+1:NOPT*I)=LCOORDS(1:NOPT)
                      IF (DEBUG) PRINT '(A,I6,A,F20.10)',' newneb> energy of image ',I,' is ',EREAL
                   ENDDO
                   DEALLOCATE(DELTAX,VNEW,LCOORDS)
                   KNOWE=.FALSE.; KNOWG=.FALSE.; KNOWH=.FALSE.
!                  REDOKADD=.FALSE.
                   REDOPATH2=.FALSE.
                ELSEIF (REDOPATH) THEN
                   XYZ(NOPT+1:NOPT*2) = TSREDO(1:NOPT)
                   PRINT '(A)','newneb> Setting tangent vector and hence initial eigenvector guess to difference between minima'
                   TANVEC(1:NOPT,1)=XYZ(1:NOPT)-XYZ(NOPT*2+1:NOPT*3)
                   EEE(2)=1.0D100
                ELSE
! IF_NOT_BIOVIA
                   IF (INTCONSTRAINTT.AND.(.NOT.(QCIDNEBT))) THEN
                      ! QCI method to construct initial band.

!                     XYZ(1:NOPT)=QQ(1:NOPT)
!                     XYZ(NOPT+1:NOPT*(NIMAGE+1))=INTNEBIMAGES(1:NOPT*NIMAGE)
                      IF (.NOT.ALLOCATED(QCIXYZ)) THEN
                         PRINT *,'newneb> ERROR *** QCIXYZ is not allocated'
                         STOP
                      ENDIF
                      IF (MOREPRINTING) CALL KEYNEBPRINT
                      XYZ(1:NOPT*(NIMAGE+2))=QCIXYZ(1:NOPT*(NIMAGE+2))
                      DEALLOCATE(QCIXYZ)
!                     PRINT '(A)','newneb> dumping qcixyz in neb.1.xyz'
!                     CALL RWG("w",.False.,1)

                      IF (RIGIDINIT) THEN
!
! Read coordinates in as atomistic xyz, but we need rigid body coordinates.
! It's possible that the xyz input file doesn't respect the rigid body constraints
! If so, converting to rigid body coordinates will restore the rigid body constraints in the
! best alignment MINPERMDIST can find to the cartesian input coords.
! This will not work very well if the 'rigid bodies' in the guessed path are very different
! to those in the reference structure.
!
                          LDEBUG = DEBUG
                          DEBUG = .FALSE. ! We can't have DEBUG set here, because it will detect that the coordinates being
                                          ! read in don't exactly match the rigid body constraints. But this is expected
                                          ! behaviour in this case.
!
! need coordinates in rigid format, but this is done below for RIGIDIINIT and not (SLERPT .OR. RIGIDMOLECULEST)
! Must NOT do it twice in a row! DJW
!
!                         PRINT '(A)','newneb> NOT calling GENRIGID_IMAGE_CTORIGID D NOT'
!                         CALL GENRIGID_IMAGE_CTORIGID(NIMAGE, XYZ) 
                          DEBUG = LDEBUG
                      ENDIF

!                     XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+2))=FINFIN(1:NOPT)
!                     CALL MAKEIMAGE(EINITIAL,EFINAL,QQ,FINFIN)
!
! Now we should respace the images to allow for the fact that the geometries of the local
! minima will lie off the interpolation between the images for the constrained potential.
! We have just put the images in XYZ so we can use INTNEBIMAGES for the respacing as in
! MAKEINTNEBIMAGE.
! We should already have the right permutational isomers for the end points, so just use
! newmindist here.
!
!                     IF (.NOT.WHOLEDNEB) THEN
                      IF (.FALSE.) THEN ! looks like we don't need any of this for QCI/LRB
                         IF (ALLOCATED(INTNEBIMAGES)) DEALLOCATE(INTNEBIMAGES)
                         ALLOCATE(INTNEBIMAGES(NIMAGE*NOPT))
                         PERMDISTSAVE=PERMDIST
                         PERMDIST=.FALSE.
                         LDEBUG = DEBUG
                         DEBUG = .FALSE. 
                          PRINT '(A)','newneb> before calling GENRIGID_IMAGE_RIGIDTOC(NIMAGE, XYZ) cant be here'
                         CALL GENRIGID_IMAGE_RIGIDTOC(NIMAGE, XYZ)
                         DEBUG = LDEBUG

                         WRITE(*,*) 'newneb> assuming QCI images are already aligned - XYZ is now in Cartesians'

!                        WRITE(*,*) "sn402: Changed NEWNEB so that it calls ALIGN_DECIDE instead of MINPERMDIST"
!                        WRITE(*,*) "At the time of writing, this WRITE statement is inside an IF(FALSE) block, ",&
!                                   "so it hasn't been tested"
!                        WRITE(*,*) "If you are reading this, please check carefully that this part of the code is working as"
!                        WRITE(*,*) "you expect, then remove this message!"
!                        CALL ALIGN_DECIDE(XYZ(NOPT+1:2*NOPT),XYZ(1:3*NATOMS),NATOMS,DEBUG, &
! &                       PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
!                        CALL ALIGN_DECIDE(XYZ(NOPT*NIMAGE+1:NOPT*(NIMAGE+1)), &
! &                       XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+2)),NATOMS,DEBUG, &
! &                       PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
                         PERMDIST=PERMDISTSAVE

                         TOTALDIST=0.0D0
                         DO J2=1,NIMAGE+1
                            XDUMMY=0.0D0
                            DO J5=1,3*NATOMS
                               XDUMMY=XDUMMY+( XYZ((J2-1)*3*NATOMS+J5) - XYZ(J2*3*NATOMS+J5) )**2
                            ENDDO
                            XDUMMY=SQRT(XDUMMY)
                            TOTALDIST=TOTALDIST+XDUMMY
                         ENDDO

!                         LDTOTAL=0.0D0
!                         DINCREMENT=0.01D0
!                         NDONE=1
!                         imageloop1: DO J2=1,NIMAGE+1
!                            XDUMMY=0.0D0
!                            DO J5=1,3*NATOMS
!                               XDUMMY=XDUMMY+( XYZ((J2-1)*3*NATOMS+J5) - XYZ(J2*3*NATOMS+J5) )**2
!                            ENDDO
!                            XDUMMY=SQRT(XDUMMY)
!                            LDIST=0.0D0
!                            DO WHILE (LDIST.LE.XDUMMY)
!                               LDIST=LDIST+DINCREMENT
!                               IF (LDIST+LDTOTAL.GE.NDONE*TOTALDIST/(NIMAGE+1)) THEN
!                                  INTNEBIMAGES(NOPT*(NDONE-1)+1:NOPT*NDONE)=((XDUMMY-LDIST)*XYZ((J2-1)*3*NATOMS+1:J2*3*NATOMS)+ &
!  &                                                          LDIST*XYZ(J2*3*NATOMS+1:(J2+1)*3*NATOMS))/XDUMMY
!                                  NDONE=NDONE+1
!                                  IF (NDONE.GT.NIMAGE) EXIT imageloop1
!                               ENDIF
!                            ENDDO
!                            LDTOTAL=LDTOTAL+XDUMMY
!                         ENDDO imageloop1
!                         XYZ(NOPT+1:NOPT*(NIMAGE+1))=INTNEBIMAGES(1:NOPT*NIMAGE)
                          DEALLOCATE(INTNEBIMAGES)
                         LDEBUG=DEBUG
                         DEBUG=.FALSE. 
                         PRINT '(A)','newneb> NOT calling GENRIGID_IMAGE_CTORIGID F'
!                        CALL GENRIGID_IMAGE_CTORIGID(NIMAGE, XYZ)
                         DEBUG=LDEBUG
                      ENDIF
                   ELSEIF (INTINTERPT) THEN
                      CALL INTINTERPOLATE(QQ,FINFIN,NINTIM,NIMAGE,X,DESMDEBUG,FAILED)
! ELSE_IF_BIOVIA
! DOBIOVIA                   IF (.FALSE.) THEN
! END_IF_NOT_BIOVIA
                   ELSE
!
! This is where we usually end up! Construct a band completely from scratch using
! (spherical) linear interpolation.
!

! IF_NOT_BIOVIA
                      IF(RIGIDINIT .AND. (RIGIDMOLECULEST .OR. SLERPT)) THEN
!
! MAKEIMAGE doesn't actually use the initial/final coordinates passed into it unless
! we have MORPHT set (it just uses the endpoints in XYZ instead) so the following
! line is usually unnecessary.
!
                        CALL MAKEIMAGE(EINITIAL,EFINAL,RIGIDQ,RIGIDFIN)
! ELSE_IF_BIOVIA
! DOBIOVIA                      IF (.FALSE.) THEN
! END_IF_NOT_BIOVIA
                      ELSE
                        CALL MAKEIMAGE(EINITIAL,EFINAL,QQ,FINFIN)
                      ENDIF
                   ENDIF
                ENDIF
             ENDIF

! hk286
! Currently, the default with RIGIDINIT is still to do unrigidified Cartesian interpolation, which allows
! the rigid bodies to deform.
! Then when we go through this block and convert to rigid body coordinates the atoms get forced
! back into their correct rigid body structures.
! If the keywords RIGIDMOLECULES (intended for a system entirely composed of identical rigid bodies)
! or SLERP are set, then instead we convert to rigid body coordinates before interpolation and do the interpolation
! using the iSLERP procedure. In that case, the following block is unnecessary.
! In the future, we may change the default behaviour so that rigid bodies are always interpolated correctly.
!
! We are in Cartesian coordinates for INTCONSTRAINTT. Need to be in rigid body coordinates for DNEB call.
!
! IF_NOT_BIOVIA
             IF (RIGIDINIT .AND. .NOT. (SLERPT .OR. RIGIDMOLECULEST)) THEN
                IF (DEBUG) PRINT '(A)','newneb> calling GENRIGID_IMAGE_CTORIGID'
                CALL GENRIGID_IMAGE_CTORIGID(NIMAGE, XYZ)
                ATOMRIGIDCOORDT = .FALSE.  !  (TRUE for Caresians, FALSE for rigid)
                IF (DEBUG) PRINT '(A,L5)','newneb> ATOMRIGIDCOORDT (TRUE for Caresians, FALSE for rigid) set to ',ATOMRIGIDCOORDT  
             ENDIF
! hk286

! preoptimise if requested        
             IF (SQVVGUESS) THEN
                STEPTOT = 5.0D0 ! TO AVOID GRADIENT SCALING DURING SQVV
                CALL NEBSQVV(NOPT*NIMAGE)
             ENDIF
! END_IF_NOT_BIOVIA

!!!!!!!!!!!!! Perform the actual NEB job here !!!!!!!!!!!!!

             NPERSIST=0
! use the c++ neb routines
! IF_NOT_BIOVIA
             IF (CPPNEBT) THEN
                IF(RIGIDINIT) THEN
! Seems strange to call this with the number of coords to
! optimise set to NOPT (which is the atomistic value)
! rather than DEGFREEDOMS. But in fact this is used to
! initialise arrays, so if we use DEGFREEDOMS then we end up
! running over the ends of arrays.
!                   PRINT '(A,L5)','newneb> calling neb_example for CPPNEBT rigidinit ATOMRIGIDCOORDT=',ATOMRIGIDCOORDT 
                    call neb_example(NIMAGE,NOPT,XYZ,EEE,NITERMAX)
                ELSE
                    CALL NEB_EXAMPLE(NIMAGE,NOPT,XYZ,EEE,NITERMAX)
                ENDIF
! ELSE_IF_BIOVIA
! DOBIOVIA                IF (.FALSE.) THEN
! END_IF_NOT_BIOVIA
             ELSEIF ((.NOT.REDOPATH).OR.REDOPATHNEB) THEN
!               PRINT '(A,L5)','newneb> should be here, ATOMRIGIDCOORDT=',ATOMRIGIDCOORDT
                SELECT CASE(MINTYPE)
                CASE("lbfgs")
! IF_NOT_BIOVIA
                   IF (UNRST) THEN
                      CALL NEBBFGSINT(NINTS*NIMAGE,NEBMUPDATE)
! ELSE_IF_BIOVIA
! DOBIOVIA                   IF (.FALSE.) THEN
! END_IF_NOT_BIOVIA
                   ELSE
!                     IF (DEBUG) PRINT '(A,L5)','newneb> calling NEBBFGS ATOMRIGIDCOORDT (T for Cartesians, F for rigid)=', &
! &                                                                      ATOMRIGIDCOORDT    
                      CALL NEBBFGS(NOPT*NIMAGE,NEBMUPDATE,NPERSIST,PERSISTENT)
                   END IF
! IF_NOT_BIOVIA
                CASE("sqvv")
                   CALL NEBSQVV(NOPT*NIMAGE)
! END_IF_NOT_BIOVIA
                END SELECT
             ENDIF

!!!!!!!!!!!!! End of main NEB job !!!!!!!!!!!!!

          ENDIF
!
! save final NEB coordinates and energy profile
!
          IF (DEBUG) THEN
             WRITE(*,'(A,F12.4)') ' newneb> mean image separation is ',SEPARATION/(NIMAGE+1)
!            DO J1=1,NIMAGE+1
!               PRINT '(A,F12.4,A,I8,A,F12.4)',' newneb> NEB k is ',NEWNEBK(J1),' for gap ',J1,' value=',DVEC(J1)
!            ENDDO
          ENDIF

          IF (DUMPNEBEOS) CALL WRITEPROFILE(0)
          IF (DUMPNEBXYZ) CALL RWG("w",.False.,0)
          IF (DUMPNEBPTS) CALL SAVEBANDCOORD

          IF (OLDCONNECT) THEN ! FIND THE HIGHEST ENERGY IMAGE
! IF_NOT_BIOVIA
             IF (DESMINT) THEN
                print*, 'newneb>> ERROR! OLDCONNECT not implemented with DESMINT'
                STOP
! ELSE_IF_BIOVIA
! DOBIOVIA             IF (.FALSE.) THEN
! DOBIOVIA
! END_IF_NOT_BIOVIA
             ENDIF
               EMAX=MAXVAL(EIMAGE)
               JMAX = 0
               DO J1=1,NIMAGE
                    IF (EMAX == EIMAGE(J1)) JMAX=J1
               ENDDO
               QQ = X(NOPT*(JMAX-1)+1:NOPT*JMAX)
          ENDIF
          
          CALL MYCPU_TIME(ENDTIME,.FALSE.)

          WRITE(*,*) "Time to go through NEB: ", ENDTIME-STARTTIME
          NMINFOUND=0
          IF (TSRESET) NTSFOUND=0
!
!  Current structure precludes searching the NEB profile for
!  both minima and ts. Could perhaps change this. Current philosophy is
!  that if we have persistent minima we should start the whole DNEB again
!  without ts searches.
!
          IF (NPERSIST.GT.0) THEN
             PTEST=.FALSE.
             IF (DEBUG) PTEST=.TRUE.
             PRINT '(A,I8)',' newneb> number of persistent minima in DNEB profile=',NPERSIST 
             DO J1=2,NIMAGE+1
                IF (PERSISTENT(J1)) THEN
                   KNOWG=.FALSE.
                   KNOWE=.FALSE. ! could use EEE value
                   IF (.NOT.ALLOCATED(VNEW)) ALLOCATE(VNEW(NOPT))
                   PRINT '(A,I8,A,F20.10)',' newneb> minimising image ',J1,' initial energy=',EEE(J1)
                   CALL MYLBFGS(NOPT,MUPDATE,XYZ(NOPT*(J1-1)+1:NOPT*J1),.FALSE., &
   &                            MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,.TRUE.,ITDONE,PTEST,VNEW,.TRUE.,.FALSE.)
                   IF (MFLAG) THEN
                      NMINFOUND=NMINFOUND+1
!
!  We have to communicate the minima found back to tryconnect using the data structure
!  set up for new transition states. 
!  Added new variable MINFOUND to allow for this check in tryconnect.
!  It seems impossible to make newneb see isnewmin and addnewmin for some reason.
!
                      ALLOCATE(MINFOUND(NMINFOUND)%COORD(NOPT))
                      MINFOUND(NMINFOUND)%COORD(1:NOPT)=XYZ(NOPT*(J1-1)+1:NOPT*J1)
                      MINFOUND(NMINFOUND)%E=EREAL
                      WRITE(987,'(I6)') NATOMS
                      WRITE(987,'(A,I5)') 'image ',J1
                      WRITE(987,'(A,3G20.10)') (ZSYM(J2), MINFOUND(NMINFOUND)%COORD(3*(J2-1)+1:3*(J2-1)+3),J2=1,NATOMS)
                   ENDIF
                   DEALLOCATE(VNEW)
                ENDIF
             ENDDO
          ELSE
             CALL PRINTSUMMARY
             CALL TSLOCATOR(TSRESET)  ! Perform the transition-state search
          ENDIF

! hk286
! IF_NOT_BIOVIA
          IF (RIGIDINIT) THEN
!
! sn402: We want to stay in RB coordinates for the moment. But for some reason we set ATOMRIGIDCOORDT
! to .TRUE. even though it isn't. I will investigate why we do this.
!                CALL GENRIGID_IMAGE_RIGIDTOC(NIMAGE, XYZ)   ! hk286 commented this line.
             ATOMRIGIDCOORDT = .TRUE.
          ENDIF
! END_IF_NOT_BIOVIA
! hk286

          NULLIFY(X,EIMAGE)
!         IF (ALLOCATED(DVEC)) DEALLOCATE(DVEC)
          IF (ALLOCATED(NEWNEBK)) DEALLOCATE(NEWNEBK)
!         IF (ALLOCATED(XYZ)) DEALLOCATE(XYZ)
!         IF (ALLOCATED(EEE)) DEALLOCATE(EEE)
!         IF (ALLOCATED(GGG)) DEALLOCATE(GGG)
          IF (ALLOCATED(TRUEGRAD)) DEALLOCATE(TRUEGRAD)
!         IF (ALLOCATED(SSS)) DEALLOCATE(SSS)
!         IF (ALLOCATED(RRR)) DEALLOCATE(RRR)
!         IF (ALLOCATED(DEVIATION)) DEALLOCATE(DEVIATION)
!         IF (ALLOCATED(TANVEC)) DEALLOCATE(TANVEC)
!         IF (ALLOCATED(STEPIMAGE)) DEALLOCATE(STEPIMAGE)
          IF (ALLOCATED(ORDERI)) DEALLOCATE(ORDERI)
          IF (ALLOCATED(ORDERJ)) DEALLOCATE(ORDERJ)
          IF (ALLOCATED(EPSALPHA)) DEALLOCATE(EPSALPHA)
          IF (ALLOCATED(DISTREF)) DEALLOCATE(DISTREF)
          IF (ALLOCATED(REPPOW)) DEALLOCATE(REPPOW)
!         DEALLOCATE(DVEC,NEWNEBK,XYZ,EEE,GGG,TRUEGRAD,SSS,RRR,DEVIATION,TANVEC,STEPIMAGE,ORDERI,ORDERJ,EPSALPHA,DISTREF,REPPOW)
          DEALLOCATE(XYZ, GGG, SSS, EEE, RRR, TANVEC, DVEC, DEVIATION, STEPIMAGE)
          DEALLOCATE(BADIMAGE,BADPEPTIDE)
          IF (FREEZENODEST) DEALLOCATE(IMGFREEZE)
! IF_NOT_BIOVIA
          IF (DESMINT) THEN
             DEALLOCATE(XCART)
             DEALLOCATE(XYZCART,GGGCART,DIHINFO)
          ENDIF
! END_IF_NOT_BIOVIA

     END SUBROUTINE NEWNEB
END MODULE NEWNEBMODULE
