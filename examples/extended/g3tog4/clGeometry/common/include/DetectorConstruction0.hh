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
/// \file common/include/DetectorConstruction0.hh
/// \brief Definition of the Common::DetectorConstruction0 class

#ifndef DetectorConstruction0_h
#define DetectorConstruction0_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "CLHEP/Units/SystemOfUnits.h"

class G4LogicalVolume;
class G4Material;
class G4GenericMessenger;

/// Simple detector construction with only a world volume

namespace Common
{

class DetectorConstruction0 : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction0(
       const G4String& materialName = "G4_AIR",
       G4double hx = 50*CLHEP::cm,
       G4double hy = 50*CLHEP::cm,
       G4double hz = 50*CLHEP::cm);
    ~DetectorConstruction0() override;

  public:
    // methods from base class
    G4VPhysicalVolume* Construct() override;

    // set methods
    void  SetMaterial(const G4String& materialName);
    void  SetDimensions(G4ThreeVector dimensions);

  private:
    void DefineCommands();

    G4GenericMessenger*  fMessenger = nullptr;
    G4String             fMaterialName;
    G4ThreeVector        fDimensions;
    G4LogicalVolume*     fWorldVolume = nullptr;
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

