      PROGRAM RESTRICTED_CORRELATION
      IMPLICIT REAL(8)(A-H,O-Z), INTEGER(4)(I-N)
C 
C     A program that generates length-dependent molar concentrations of
C     Worm-Like Chain (WLC) ends at four points PA, PB, PC, PD, for WLCs
C     constrained to Z>0, starting at X=0, Y=0, Z=0.
C 
C     Depends on the procedures "WLC_simulation_posZ_stepwise.f" and
C     "inv_Langevin_function.f"
C 
C     ------------------------------------------------------------------
C 
C         Edit the values of the following 5+12 parameters
C         before using the command "make concentrations"
C 
      PARAMETER(PL=0.50D0)      ! Persistence length (nm)
      PARAMETER(BL=0.36D0)      ! Contour length per monomer (nm)
      PARAMETER(NM=260)         ! Number of monomers in the polymer
      PARAMETER(NWLC=100000000) ! Number of WLC simulated; greater
C                               !    NWLS produces smoother distrubution 
      PARAMETER(TOROIDR=2.0D-1) ! Cutoff distance from the points PA, PB,
C                               ! PC, PD; equivalent of bin width in nm.
C 
C         Coordinates X, Y, Z for the points PA, PB, PC, PD:
C
      PARAMETER(PAX= 2.08612D0, PAY=-1.14450D0, PAZ= 3.59757D0)
      PARAMETER(PBX= 4.00242D0, PBY=-0.62052D0, PBZ= 1.99452D0)
      PARAMETER(PCX= 5.74549D0, PCY=-2.73123D0, PCZ= 1.40088D0)
      PARAMETER(PDX= 3.75826D0, PDY=-6.71115D0, PDZ= 2.01347D0)
C
C         End of the parameters that can be edited.
C 
      PARAMETER(TWOPI=6.283185307179586232D0) ! Do not change this parameter
      PARAMETER(AVOGADRO=6.02214076D23)       !mol^-1
      PARAMETER(UNITVOLUME=1.0D-24)           !liters per nm^3
      REAL(8) milliMolar
      DIMENSION IBINA(NM), IBINB(NM), IBINC(NM), IBIND(NM)
C 
      PARAMETER(ISEEDDIM=8)     ! Do not change this parameter
      DIMENSION ISEED(ISEEDDIM)
C 
C     The eight integers below represent the seed for the random number
C     generator; change these numbers if you want a different set of WLCs
C     to be generated. If NWLC is large enough, nothing will chamne when
C     you change the eight integers below.
C 
      DATA ISEED /293753205, 945279830, 847569340, 752095743,
     &            206706445, 362953259, 572434685, 856872232/
C 
      CALL RANDOM_SEED(SIZE=JSEEDDIM)
      IF(JSEEDDIM.GT.ISEEDDIM)THEN
        WRITE(*,*)JSEEDDIM
        STOP 'ERROR: ISEEDDIM IS TOO SMALL.'
      ENDIF
      CALL RANDOM_SEED(PUT=ISEED)
C 
      milliMolar=1.0D3/(UNITVOLUME*AVOGADRO)
      TOROIDR2=TOROIDR*TOROIDR
      PAH=SQRT(PAX*PAX+PAY*PAY)
      PBH=SQRT(PBX*PBX+PBY*PBY)
      PCH=SQRT(PCX*PCX+PCY*PCY)
      PDH=SQRT(PDX*PDX+PDY*PDY)
C 
      DO 10,I=1,NM
      IBINA(I)=0
      IBINB(I)=0
      IBINC(I)=0
      IBIND(I)=0
  10  CONTINUE
C 
      DO 100,IWLC=1,NWLC
      CALL WLC_SIMULATION_POSZ_SETUP(X,Y,Z,PL,BL)
      DO 110,I=1,NM
      CALL WLC_SIMULATION_POSZ_ONESTEP(X,Y,Z)
      H=SQRT(X*X+Y*Y)
      DIFFH=H-PAH
      DIFFZ=Z-PAZ
      IF(DIFFH*DIFFH+DIFFZ*DIFFZ.LE.TOROIDR2)IBINA(I)=IBINA(I)+1
      DIFFH=H-PBH
      DIFFZ=Z-PBZ
      IF(DIFFH*DIFFH+DIFFZ*DIFFZ.LE.TOROIDR2)IBINB(I)=IBINB(I)+1
      DIFFH=H-PCH
      DIFFZ=Z-PCZ
      IF(DIFFH*DIFFH+DIFFZ*DIFFZ.LE.TOROIDR2)IBINC(I)=IBINC(I)+1
      DIFFH=H-PDH
      DIFFZ=Z-PDZ
      IF(DIFFH*DIFFH+DIFFZ*DIFFZ.LE.TOROIDR2)IBIND(I)=IBIND(I)+1
 110  CONTINUE
c 
      kwlc=nwlc/1000
      jwlc=iwlc/kwlc
      if(jwlc*kwlc.eq.iwlc)write(*,'(F5.1,"%")')0.1D0*dble(jwlc)
c 
 100  CONTINUE
C 
      TOROIDVOLUME=(0.5D0*TWOPI*TOROIDR2)*(TWOPI*PAH)
      OPEN(1,FILE='concentration_PA.txt')
      DO 210,I=1,NM
      Y=DBLE(IBINA(I))/(TOROIDVOLUME*DBLE(NWLC))
      WRITE(1,*)I,Y*milliMolar
 210  CONTINUE
      CLOSE(1)
C 
      TOROIDVOLUME=(0.5D0*TWOPI*TOROIDR2)*(TWOPI*PBH)
      OPEN(1,FILE='concentration_PB.txt')
      DO 220,I=1,NM
      Y=DBLE(IBINB(I))/(TOROIDVOLUME*DBLE(NWLC))
      WRITE(1,*)I,Y*milliMolar
 220  CONTINUE
      CLOSE(1)
C 
      TOROIDVOLUME=(0.5D0*TWOPI*TOROIDR2)*(TWOPI*PCH)
      OPEN(1,FILE='concentration_PC.txt')
      DO 230,I=1,NM
      Y=DBLE(IBINC(I))/(TOROIDVOLUME*DBLE(NWLC))
      WRITE(1,*)I,Y*milliMolar
 230  CONTINUE
      CLOSE(1)
C 
      TOROIDVOLUME=(0.5D0*TWOPI*TOROIDR2)*(TWOPI*PDH)
      OPEN(1,FILE='concentration_PD.txt')
      DO 240,I=1,NM
      Y=DBLE(IBIND(I))/(TOROIDVOLUME*DBLE(NWLC))
      WRITE(1,*)I,Y*milliMolar
 240  CONTINUE
      CLOSE(1)
C 
      STOP
      END

