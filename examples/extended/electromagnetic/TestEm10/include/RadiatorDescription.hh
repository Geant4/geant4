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
/// \file electromagnetic/TestEm10/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
//
//

#ifndef RadiatorDescription_h
#define RadiatorDescription_h 1

#include "globals.hh"

class G4LogicalVolume;
class G4Material;

struct RadiatorDescription
{
  RadiatorDescription()
    : fLogicalVolume(0), fFoilMaterial(0), fGasMaterial(0), 
      fFoilThickness(0.), fGasThickness(0.), fFoilNumber(0) {}

  G4LogicalVolume*  fLogicalVolume;
  G4Material*       fFoilMaterial;
  G4Material*       fGasMaterial;
  G4double          fFoilThickness;
  G4double          fGasThickness;
  G4int             fFoilNumber;
};

#endif


