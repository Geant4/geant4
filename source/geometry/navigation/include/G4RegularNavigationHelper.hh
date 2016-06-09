//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4RegularNavigationHelper.hh,v 1.1 2009-01-27 09:31:29 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4RegularNavigationHelper
//
// Class description:
//
// Utility class for navigation on regular structures, providing step
// lengths counting for each regular voxel of the structure.

// Author: Pedro Arce, November 2008
// --------------------------------------------------------------------
#ifndef G4RegularNavigationHelper_HH
#define G4RegularNavigationHelper_HH

#include <vector>
#include "globals.hh"

class G4RegularNavigationHelper
{
  public:

    G4RegularNavigationHelper();
   ~G4RegularNavigationHelper();
  
    static void ClearStepLengths();
    static void AddStepLength( G4int copyNo, G4double slen );

    static std::vector< std::pair<G4int,G4double> > theStepLengths;
};

#endif
