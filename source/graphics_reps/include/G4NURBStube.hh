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
// $Id$
//
// 
// Olivier Crumeyrolle  12 September 1996

// Tube builder prototype
// OC 090796

// Class Description:
// Tube builder prototype for NURBS.
// See documentation in graphics_reps/doc for details.

#ifndef __C_G4NURBStube__
#define __C_G4NURBStube__ 1 

#include "G4NURBS.hh"

class  G4NURBStube : public G4NURBS
{
  public:
    G4NURBStube(G4double RMIN, G4double RMAX, G4double DZ);
    const char*  Whoami() const;
};
#endif
// end of __C_G4NURBStube__
