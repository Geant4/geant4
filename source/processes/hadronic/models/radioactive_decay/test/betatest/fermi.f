*    
*     ********************************************************************
*     * License and Disclaimer                                           *
*     *                                                                  *
*     * The  Geant4 software  is  copyright of the Copyright Holders  of *
*     * the Geant4 Collaboration.  It is provided  under  the terms  and *
*     * conditions of the Geant4 Software License,  included in the file *
*     * LICENSE and available at  http://cern.ch/geant4/license .  These *
*     * include a list of copyright holders.                             *
*     *                                                                  *
*     * Neither the authors of this software system, nor their employing *
*     * institutes,nor the agencies providing financial support for this *
*     * work  make  any representation or  warranty, express or implied, *
*     * regarding  this  software system or assume any liability for its *
*     * use.  Please see the license in the file  LICENSE  and URL above *
*     * for the full disclaimer and the limitation of liability.         *
*     *                                                                  *
*     * This  code  implementation is the result of  the  scientific and *
*     * technical work of the GEANT4 collaboration.                      *
*     * By using,  copying,  modifying or  distributing the software (or *
*     * any work based  on the software)  you  agree  to acknowledge its *
*     * use  in  resulting  scientific  publications,  and indicate your *
*     * acceptance of all terms of the Geant4 Software license.          *
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
