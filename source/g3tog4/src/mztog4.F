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
*    
*     $Id: mztog4.F 67982 2013-03-13 10:36:03Z gcosmo $
*    
      subroutine mztog4
************************************************************************
*
*       mztog4
*
*       Within a Geant3 program, call this routine to perform conversion
*       to Geant4 on G3 data structures in memory.
*
*...  History:
*     15-Jan-1997 Lockman : Put tog4 in a separate file
************************************************************************
      implicit none
#include "G3toG4.inc"
*
***     conversion initialization
      dogeom = .false.
      context = '----'
      lunlist = 98
      open(unit=lunlist,file='g3calls.dat',status='unknown')
      luncode = 99
      nfile = 1
      call g3source
*
***     do conversion
      call tog4
*
      end
*
