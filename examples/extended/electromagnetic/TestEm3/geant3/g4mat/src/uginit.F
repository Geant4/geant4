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

      SUBROUTINE UGINIT                                                         
*                                                                               
*     To initialise GEANT/USER  program and read data cards                     
*
#include "geant321/gcphys.inc"
#include "geant321/gccuts.inc"
#include "ugmate.inc"
#include "ggvalues.inc"
*
      CHARACTER*20 filnam                                                                                                                                                                                                      
*                                                                               
* *** Define the GEANT parameters
      CALL GINIT 
                                                                          
* *** read data cards
      PRINT *, 'G3 > gives the filename of the data cards to be read:'
      READ (*,'(A)') filnam
      IF (filnam.EQ.' ') filnam = 'run01.dat'
      OPEN (unit=5,file=filnam,status='unknown',form='formatted')
*
* *** define data cards
      call FFKEY ('MATE' ,Imat,1,'INTEGER')
      call FFKEY ('G4VAL',g4value(1),200,'REAL')  
*                                                                    
* *** read data cards                                                           
      CALL GFFGO 
*
* *** overwrite some initialization
      iloss  = 4
      ppcutm = cutele
*                                                                                      
      CALL GZINIT                                                               
      CALL GPART
*                                                                               
* *** materials definition
      CALL UGIMAT                                                               
*                                                                               
* *** Energy loss and cross-sections initialisations
      CALL GPHYSI                                                               
*                                                                               
      CALL GPRINT('MATE',0)
*
* *** ready for geant3/geant4 comparison
      call UG3G4      
*                                                                               
      END                                                                       
