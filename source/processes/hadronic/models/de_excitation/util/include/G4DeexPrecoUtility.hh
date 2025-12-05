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
// Geant4 header G4DeexPrecoUtility
//
// Author V.Ivanchenko 19.05.2025
//
// Utilities used at initialisation of the de-excitation module
//

#ifndef G4DeexPrecoUtility_h
#define G4DeexPrecoUtility_h 1

#include "globals.hh"

class G4DeexPrecoUtility
{
public:

  // Compute correction factor
  // Data comes from Dostrovsky, Fraenkel and Friedlander
  // Physical Review, vol 116, num. 3 1959

  static G4double CorrectionFactor(const G4int index, const G4int Z,
				   const G4double A13,
				   const G4double bCoulomb,
				   const G4double ekin);

  static G4double ProtonKValue(const G4int Z);

  static G4double AlphaKValue(const G4int Z);

  static G4double ProtonCValue(const G4int Z);

  static G4double AlphaCValue(const G4int Z);

  static G4double LevelDensity(const G4int Z, const G4int A, const G4int index);

};

#endif


