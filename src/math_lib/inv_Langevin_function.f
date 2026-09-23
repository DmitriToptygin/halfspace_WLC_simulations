      FUNCTION INV_LANGEVIN(X)
      REAL(8) INV_LANGEVIN, X
      INTENT(IN) X
C 
C     Inverse Langevin function: INV_LANGEVIN(x)
C 
C     Langevin function:  LANGEVIN(y)=coth(y)-1/y
C 
C     if x=LANGEVIN(y), then y=INV_LANGEVIN(x)
C 
C     This procedure uses Newton-Raphson method to calculate the inverse
C     of the Langevin function. Application of the Newton-Raphson method
C     requires calculating Langevin function and its derivative.
C 
C     Langevin function suffers from catastrophic cancellation near y=0
C     because both coth(y) and 1/y tend to infinity as y approaches 0,
C     however, the difference between coth(y) and 1/y tends to 0 as y
C     approaches 0. Thus, subtracting 1/y from coth(y) results in poor
C     accuracy near y=0. To attain full REAL(8) accuracy, this procedure
C     uses a continued fraction expansion of Langevin function for y
C     values within the range from -3 to +3. The continued fraction
C     expansion of Langevin function has been derived from the Lambert's
C     continued fraction expansion of tanh(y). The derivative of
C     Langevin function also suffers from catastrophic cancellation near
C     y=0. This procedure utilizes an analytical derivative of the
C     continued fraction expansion of Langevin function for y values
C     within the range from -3 to +3.
C 
C     ------------------------------------------------------------------
C     Copyright (c) 2026 Dmitri Toptygin
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
      REAL(8) X2, Y, TY, Y2, E, CH, LANGEVIN, DERLANGEVIN, D,P,Q, DELTA
      X2=X*X
      IF(X2.GE.1.0D0)THEN            ! checking whether X is within the domain.
        Y=0.0D0
        INV_LANGEVIN=Y/Y             ! If X is outside the domain, then
        GOTO 20                      !  the function returns a NaN value
      ENDIF
      Y=X*(3.0D0-X2)/(1.0D0-X2)      ! starting approximation by A. Cohen (1991).
C 
  10  CONTINUE                       ! Newton-Raphson method starts here
C 
C     Accurate calculation of LANGEVIN(y)
C                                               d LANGEVIN(y)
C                         and DERLANGEVIN(y) = ---------------
C                                                   d y
      TY=Y+Y
      Y2=Y*Y
      IF    (Y.LT.-3.0D0)THEN
        E=EXP( TY)
        CH=(E+1.0D0)/(E-1.0D0)
        LANGEVIN=CH-1.0D0/Y
        DERLANGEVIN=(1.0D0-CH)*(1.0D0+CH)+1.0D0/Y2
      ELSEIF(Y.LE. 3.0D0)THEN
        D=0.21727302006229138D0
        P=25.329582258975275D0
        Q=Y2/P
        D=(TY-D*Q)/P
        P=Q+23.0D0
        Q=Y2/P
        D=(TY-D*Q)/P
        P=Q+21.0D0
        Q=Y2/P
        D=(TY-D*Q)/P
        P=Q+19.0D0
        Q=Y2/P
        D=(TY-D*Q)/P
        P=Q+17.0D0
        Q=Y2/P
        D=(TY-D*Q)/P
        P=Q+15.0D0
        Q=Y2/P
        D=(TY-D*Q)/P
        P=Q+13.0D0
        Q=Y2/P
        D=(TY-D*Q)/P
        P=Q+11.0D0
        Q=Y2/P
        D=(TY-D*Q)/P
        P=Q+ 9.0D0
        Q=Y2/P
        D=(TY-D*Q)/P
        P=Q+ 7.0D0
        Q=Y2/P
        D=(TY-D*Q)/P
        P=Q+ 5.0D0
        Q=Y2/P
        D=(TY-D*Q)/P
        P=Q+ 3.0D0
        Q=Y/P
        D=(1.0D0-D*Q)/P
        LANGEVIN=Q
        DERLANGEVIN=D
      ELSE
        E=EXP(-TY)
        CH=(1.0D0+E)/(1.0D0-E)
        LANGEVIN=CH-1.0D0/Y
        DERLANGEVIN=(1.0D0-CH)*(1.0D0+CH)+1.0D0/Y2
      ENDIF
C 
C     Accurate calculation of LANGEVIN(y) and DERLANGEVIN(y) ends here,
C     Newton-Raphson method continues
C 
      DELTA=X-LANGEVIN
      Y=Y+DELTA/DERLANGEVIN
      IF(ABS(DELTA).GT.1.0D-13)GOTO 10
      INV_LANGEVIN=Y
  20  RETURN
      END

