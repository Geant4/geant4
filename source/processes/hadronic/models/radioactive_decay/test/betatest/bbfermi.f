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
      FUNCTION BBFERMI(X)
      Z0 = 100
      PI=4.0*ATAN(1.0)
      P=SQRT(X*X-1.0) 
      U=Z0/137.0
      S=SQRT(1.0-U*U) - 1.
      Y = 2*PI*U*SQRT(1+P*P)/P
      A1 = U*U*X*X + P*P/4.
      A2 = abs(Y/(1-exp(-Y)))
      BBFERMI=A1**S * A2
      print *, A1,A2,BBFERMI
      RETURN
      END
