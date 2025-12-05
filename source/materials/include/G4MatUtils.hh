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
// Geant4 header G4MatUtils
//
// Author V.Ivanchenko 30.09.2025
//
// Utilities used at initialisation of material properties
//

#ifndef G4MatUtils_HH
#define G4MatUtils_HH

#include "globals.hh"
#include "G4ExtendedPhysicsVector.hh"

class G4MatUtils
{
public:

  // Build G4ExtendedPhysicsVector using data provided in the fixed format
  // of length = number of strings with the same data:  
  //    energy and N partial cross sections
  static G4ExtendedPhysicsVector*
  BuildExtendedVector(const G4String& dirpath, const G4String& filename,
		      const G4int N, const G4int length,
		      const G4double unitE = 1.0, const G4double unitS = 1.0);

};

#endif
