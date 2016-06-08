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

      SUBROUTINE UG3G4                                                         
*
*                                                                               
#include "geant321/gcmulo.inc"
#include "geant321/gckine.inc"
#include "ugmate.inc"
#include "ggvalues.inc"
*
*
* *** fetch Geant3 tables      
      call gftmat(imat,ikine,'LOSS',nekbin,elow,g3value,dum,idum)
*
* *** histograms to plot the difference
      elow(nek1) = ekmax
      call HBOOKB (1,'(g4-g3)/g3 in percent',nekbin,elow,0.)
      call HBOOKB (2,'abs((g4-g3)/g3) in percent',nekbin,elow,0.)           
*
* *** compute difference
      do i=1,nekbin
        dif = 100*(g4value(i)-g3value(i))/g3value(i)
	call hfill (1, elow(i), 0., dif)
	call hfill (2, elow(i), 0., abs(dif))	
        print *,i,elow(i),g4value(i),g3value(i), 'dif(%)=',dif
      enddo      
*                                                                               
      END                                                                       
