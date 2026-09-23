      PROGRAM GENERATE_DISTANCE_DISTRIBUTION_POSZ
      IMPLICIT REAL(8)(A-H,O-Z), INTEGER(4)(I-N)
C 
C     A program that generates end-to-end distance distribution for a
C     Worm-Like Chain (WLC) constrained to Z>0, starting at Z=0.
C 
C     Depends on the procedures "WLC_simulation_posZ_endpoint.f" and
C     "inv_Langevin_function.f"
C 
C     ------------------------------------------------------------------
C 
C         Edit the values of the following six parameters
C         before using the command "make distance_distribution"
C 
      PARAMETER(PL=0.50D0)      ! Persistence length (nm)
      PARAMETER(BL=0.36D0)      ! Contour length per monomer (nm)
      PARAMETER(NM=259)         ! Number of monomers in the polymer
      PARAMETER(NWLC=100000000) ! Number of WLC simulated; greater
C                               !    NWLS produces smoother distrubution 
      PARAMETER(DR=0.04D0)      ! Bin width for the distance distribution
      PARAMETER(NBIN=10000)     ! Total number of bins
      DIMENSION IBIN(NBIN)
C 
      PARAMETER(ISEEDDIM=8)     ! Do not change this parameter
      DIMENSION ISEED(ISEEDDIM)
C 
C     The eight integers below represent the seed for the random number
C     generator; change these numbers if you want a different set of WLCs
C     to be generated. If NWLC is large enough, nothing will chamne when
C     you change the eight integers below.
C 
      DATA ISEED /293753259, 945279830, 847569340, 752095743,
     &            206706445, 362953259, 572434685, 856872232/
C 
      CALL RANDOM_SEED(SIZE=JSEEDDIM)
      IF(JSEEDDIM.GT.ISEEDDIM)THEN
        WRITE(*,*)JSEEDDIM
        STOP 'ERROR: ISEEDDIM IS TOO SMALL.'
      ENDIF
      CALL RANDOM_SEED(PUT=ISEED)
C 
      DO 10,N=1,NBIN
      IBIN(N)=0
  10  CONTINUE
C 
      MBIN=0
      DO 20,IWLC=1,NWLC
      CALL WLC_SIMULATION_POSZ(X,Y,Z,PL,BL,NM)
      R=SQRT(X*X+Y*Y+Z*Z)
      N=NINT(0.5D0+R/DR)
      IF(N.LE.NBIN)IBIN(N)=IBIN(N)+1
      MBIN=MAX(MBIN,N)
c 
      kwlc=nwlc/1000
      jwlc=iwlc/kwlc
      if(jwlc*kwlc.eq.iwlc)write(*,'(F5.1,"%")')0.1D0*dble(jwlc)
c 
  20  CONTINUE
      MBIN=MIN(MBIN+1,NBIN)
C 
      OPEN(1,FILE='distribution.txt')
      DO 30,N=1,MBIN
      WRITE(1,*)DR*(DBLE(N)-0.5D0),DBLE(IBIN(N))/(DR*DBLE(NWLC))
  30  CONTINUE
      CLOSE(30)
C 
      STOP
      END

