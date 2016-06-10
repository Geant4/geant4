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
*     $Id: G3toG4.inc 67982 2013-03-13 10:36:03Z gcosmo $
*    
      integer maxlines
      parameter (maxlines = 750) ! max no. of cmd lines per C++ source file

      integer irot_index
      integer imat_index
      character context*4
      logical doclist      ! Build call list
      logical docode       ! Build C++ code
      logical dogeom       ! Build the Geant geometry as well
      integer lunlist      ! Output logical unit for call list
      integer luncode      ! Output logical unit for code
      integer nfile        ! Current source code file number
      integer nlines       ! Current source code line count
      common /g3tog4_c/ irot_index, imat_index,
     +  doclist, docode, dogeom, lunlist, luncode, context, nfile,
     +  nlines
