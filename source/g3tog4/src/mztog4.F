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
*    
*     $Id: mztog4.F,v 1.3 2001-07-11 09:59:03 gunter Exp $
*     GEANT4 tag $Name: not supported by cvs2svn $
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
