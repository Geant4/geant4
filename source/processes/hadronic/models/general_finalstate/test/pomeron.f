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
         FUNCTION pomeron(X)
         
         double precision impactsquare
         
         double precision Eikonal
         double precision Lambda
         double precision Z
         double precision Power
         double precision s

         double precision hbarc_squared, GeV, fermi2
         double precision pomeron_Gamma, pomeron_S, pomeron_Alpha
         double precision pomeron_Rsquare, pomeron_Alphaprime, pomeron_C
         
         double precision fac0, fac1, fac2, fac3
         
         fermi2 = 1.0E-25
         GeV = 1000.0
         hbarc_squared = 3.8937966265481e-20
         
         impactsquare = x*x*fermi2
         s = 500*GeV*GeV
          
C         pomeron_Gamma = (2.6+3.96)/GeV/GeV
         pomeron_Gamma = 3.96/GeV/GeV
C         pomeron_S = 2.7*GeV*GeV
         pomeron_S = 3.0*GeV*GeV
C         pomeron_Alpha = 0.9808
         pomeron_Alpha = 1.0808
C         pomeron_Rsquare = 1.7*3.56/GeV/GeV
         pomeron_Rsquare = 13.56/GeV/GeV
C         pomeron_Alphaprime = 1.6*0.25/GeV/GeV
         pomeron_Alphaprime = .25/GeV/GeV
C         pomeron_C = 1.4
         pomeron_C = 1.4
         
  	 fac0 = s/pomeron_S
  	 fac1 = pomeron_Alpha -1
  	 fac2 = fac0**fac1
  	 Power = pomeron_Gamma *fac2
	 Lambda = pomeron_Rsquare+pomeron_Alphaprime*log(s/pomeron_S)
	 Z = 2*pomeron_C * Power / Lambda
	 fac0 = -impactsquare/(4*Lambda*hbarc_squared)
	 Eikonal = Z/2 * exp(fac0)
	 
	 pomeron = 2/pomeron_C*(1-exp(-1*Eikonal))
	 print *,'testhpw ', pomeron
	 	
         END
