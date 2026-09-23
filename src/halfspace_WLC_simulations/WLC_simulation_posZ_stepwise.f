      SUBROUTINE WLC_SIMULATION_POSZ_SETUP(X,Y,Z,PL,BL)
      REAL(8) X, Y, Z, PL, BL
      INTENT(IN) PL, BL
      INTENT(INOUT) X, Y, Z
C 
C     Simulate a 3D worm-like chain (WLC) polymer constrained to Z > 0.
C 
C     Inputs:
C       PL  - Persistence length (nm)
C       BL  - Contour length per monomer (nm)
C       NM  - Number of monomers in the polymer (not a formal argumant).
C 
C     Inputs/Outputs:
C       X, Y, Z - Coordinates of the polymer chain monomer (Z>0)
C 
C     This prosidure is similar to WLC_SIMULATION_POSZ(X,Y,Z,PL,BL,NM),
C     but the simulation is run one step (one monomer) at a time, which
C     allows the calling program to save the current values of X, Y, Z
C     after every step rather than just after NM steps.
C 
C     The first call must be to WLC_SIMULATION_POSZ_SETUP(X,Y,Z,PL,BL);
C     followed by NM calls to WLC_SIMULATION_POSZ_ONESTEP(X,Y,Z).
C 
C     Note: the values of X, Y, Z, must not be changed between the calls
C     to WLC_SIMULATION_POSZ_SETUP() and WLC_SIMULATION_POSZ_ONESTEP().
C 
C     WLC theory:
C 
C     WLC consists of NM segments of length BL; the directoion of each
C     segment is parallel to a unit vector D=(Dx,Dy,Dz), which evolves
C     from one segment to the next one according to the following rules:
C 
C     The angle θ is defined such that the inner product of the previous
C     vector D times the next vector D equals cos(θ).
C 
C     If θ=0, then the old and the new vector are parallel and pointing
C     in the same direction; if θ=π, then they are parallel and pointing
C     in opposite directions.
C 
C     P(θ) is the probability density of angle θ,
C 
C     P(θ) = sin(θ) * exp[u*cos(θ)] * u/[exp(u)-exp(-u)]
C 
C     The ensemble mean cos(θ) is a function of the parameter u:
C 
C     ⟨cos(θ)⟩ = L(u)
C 
C     where L(...) is the Langevin function, L(x) = coth(x) - 1/x,
C 
C     where coth(...) is a hyperbolic cotangent.
C 
C     The ensemble mean cos(θ) can be expressed in terms of the contour
C     length per monomer (bond length BL) and persistence length (PL):
C 
C     ⟨cos(θ)⟩ = exp(-BL/PL)
C 
C     This gives a method for calculating parameter u:
C 
C     u = inv_Langevin[exp(-BL/PL)]
C 
C     where inv_Langevin(...) is the inverse Langevin function
C 
C     To generate random angles θ that are distributed in accordance
C     with the above probability density P(θ) it is convenient to use
C     the cumulative distribution function for cos(θ):
C 
C     CDF[cos(θ)] = {exp(u)-exp[u*cos(θ)]}/[exp(u)-exp(-u)]
C 
C     If RND is a random variable the values of which are uniformly
C     distributed in the range from 0 to +1, then the random values of
C     cos(θ) can be calculated form RND as follows:
C 
C     cos(θ) = ln{exp(u)-RND*[exp(u)-exp(-u)]}/u
C 
C     To calculate the direction of the new vector D from that of the
C     old vector D one must know not only the angle θ, but also the
C     dihedral angle formed by three segments (the one before the
C     previous, the previous, and the new one). However, this dihedral
C     angle is uniformly distributed in the range from -π to +π; in
C     other words, the dihedral angle is completely random, which makes
C     the direction of the vector before the previous one irrelevant.
C 
C     The starting direction of vector D is (0,0,1).
C 
C     The starting point of the WLC is (0,0,0).
C 
C     ------------------------------------------------------------------
C     Copyright (c) 2026 Christian M. Kaiser and Dmitri Toptygin
C 
C     Permission is hereby granted, free of charge, to any person obtaining a
C     copy of this software and associated documentation files (the "Software"),
C     to deal in the Software without restriction, including without limitation
C     the rights to use, copy, modify, merge, publish, distribute, sublicense,
C     and/or sell copies of the Software, and to permit persons to whom the
C     Software is furnished to do so, subject to the following conditions:
C 
C     The above copyright notice and this permission notice shall be included
C     in all copies or substantial portions of the Software.
C 
C     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
C     IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
C     FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
C     THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
C     LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
C     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
C     DEALINGS IN THE SOFTWARE.
C     ------------------------------------------------------------------
      REAL(8) TWOPI, BL0, CORR_FAC, U, EXPU, DIFFEU, INV_LANGEVIN
      REAL(8) DX,DY,DZ, FX,FY,FZ, GX,GY,GZ, RX,RY,RZ, DPR, C
      REAL(8) XOLD,YOLD,ZOLD, DXOLD,DYOLD,DZOLD, DX2,DY2,DZ2
      REAL(8) RND, COSTHETA, SINTHETA, PHI, SINPHI, COSPHI
      SAVE BL0, CORR_FAC, U, EXPU, DIFFEU, DX, DY, DZ
      EXTERNAL INV_LANGEVIN
      PARAMETER(TWOPI=6.283185307179586232D0)
      BL0=BL
      CORR_FAC=EXP(-BL/PL)
      U=INV_LANGEVIN(CORR_FAC)
      EXPU=EXP(U)
      DIFFEU=EXPU-EXP(-U)
      X=0.0D0
      Y=0.0D0
      Z=0.0D0
      DX=0.0D0
      DY=0.0D0
      DZ=1.0D0
      RETURN
C     ------------------------------------------------------------------
      ENTRY WLC_SIMULATION_POSZ_ONESTEP(X,Y,Z)
C 
      XOLD=X
      YOLD=Y
      ZOLD=Z
      DXOLD=DX
      DYOLD=DY
      DZOLD=DZ
 110  CONTINUE
      CALL RANDOM_NUMBER(RND)
      COSTHETA=LOG(EXPU-RND*DIFFEU)/U
      SINTHETA=SQRT(1.0D0-COSTHETA*COSTHETA)
C 
C     Constructing a unit vector (FX,FY,FZ) normal to vector (DX,DY,DZ)
C     assuming that (DX,DY,DZ) is a unit vector:
C 
      DX2=DX*DX
      DY2=DY*DY
      DZ2=DZ*DZ
      IF    (DX2.LE.DY2.AND.DX2.LE.DZ2)THEN
        FX=DX2-1.0D0
        FY=DX*DY
        FZ=DX*DZ
      ELSEIF(DY2.LE.DX2.AND.DY2.LE.DZ2)THEN
        FX=DY*DX
        FY=DY2-1.0D0
        FZ=DY*DZ
      ELSE
        FX=DZ*DX
        FY=DZ*DY
        FZ=DZ2-1.0D0
      ENDIF
      DPR=FX*FX+FY*FY+FZ*FZ
      C=1.0D0/SQRT(DPR)
      FX=C*FX
      FY=C*FY
      FZ=C*FZ
C 
C     Constructing a unit vector (GX,GY,GZ) normal to vectors (DX,DY,DZ)
C     and to vector (FX,FY,FZ) assuming that (DX,DY,DZ) is a unit vector
C     and (FX,FY,FZ) is also a unit vector normal to vector (DX,DY,DZ):
C 
      GX=DY*FZ-DZ*FY
      GY=DZ*FX-DX*FZ
      GZ=DX*FY-DY*FX
C 
C     Generating a random unit vector (RX,RY,RZ) that is normal to
C     vector (DX,DY,DZ); all directions of (RX,RY,RZ) in the plane
C     normal to vector (DX,DY,DZ) have eqial probability:
C 
      CALL RANDOM_NUMBER(RND)
      PHI=TWOPI*(RND-0.5D0)
      SINPHI=SIN(PHI)
      COSPHI=COS(PHI)
      RX=COSPHI*FX+SINPHI*GX
      RY=COSPHI*FY+SINPHI*GY
      RZ=COSPHI*FZ+SINPHI*GZ
C 
C     Generating new vector (DX,DY,DZ) from the old one and from the
C     random vector (RX,RY,RZ)
C 
      DX=COSTHETA*DX+SINTHETA*RX
      DY=COSTHETA*DY+SINTHETA*RY
      DZ=COSTHETA*DZ+SINTHETA*RZ
C 
C     Normalizing new vector (DX,DY,DZ) to unit length to prevent
C     accumulation of numerical errors
C 
      DPR=DX*DX+DY*DY+DZ*DZ
      C=1.0D0/SQRT(DPR)
      DX=C*DX
      DY=C*DY
      DZ=C*DZ
      X=X+BL0*DX
      Y=Y+BL0*DY
      Z=Z+BL0*DZ
      IF(Z.LE.0.0D0)THEN
        X=XOLD
        Y=YOLD
        Z=ZOLD
        DX=DXOLD
        DY=DYOLD
        DZ=DZOLD
        GOTO 110
      ENDIF
C 
      RETURN
      END

