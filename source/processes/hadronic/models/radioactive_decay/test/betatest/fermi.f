*    
*     ********************************************************************
*     * DISCLAIMER                                                       *
*     *                                                                  *
*     * The following disclaimer summarizes all the specific disclaimers *
*     * of contributors to this software. The specific disclaimers,which *
*     * govern, are listed with their locations in:                      *
*     *   http://cern.ch/geant4/license                                  *
*     *                                                                  *
*     * Neither the authors of this software system, nor their employing *
*     * institutes,nor the agencies providing financial support for this *
*     * work  make  any representation or  warranty, express or implied, *
*     * regarding  this  software system or assume any liability for its *
*     * use.                                                             *
*     *                                                                  *
*     * This  code  implementation is the  intellectual property  of the *
*     * GEANT4 collaboration.                                            *
*     * By copying,  distributing  or modifying the Program (or any work *
*     * based  on  the Program)  you indicate  your  acceptance of  this *
*     * statement, and all its terms.                                    *
*     ********************************************************************
*    
      FUNCTION FERMI(X)
      Z0 = 35
      A0 = 71
      PI=4.0*ATAN(1.0)
      P=SQRT(X*X-1.0) 
      BETA=P/X
      U=Z0/137.0
      S=SQRT(1.0-U*U)
      R=1.48*EXP((1./3.)*ALOG(A0))
      RO=R/386.157
      ETA=U/BETA
      FI=ATAN(S/ETA)
      S2ETA2=S*S+ETA*ETA
      A1=4.0*PI*(1+S)
      A2=GAMMA(2.0*S)
      A2=A2*A2
      A3=EXP(2.0*(S-1.0)*ALOG(2.0*P*RO))
      A4=EXP((S-0.5)*ALOG(S2ETA2))
      A5=2.0*FI*ETA-2.0*S+S/6.0/S2ETA2
      FERMI=A1/A2*A3*A4*EXP(A5)
      print *, A1,A2,A3,A4,A5,FI
      RETURN
      END
