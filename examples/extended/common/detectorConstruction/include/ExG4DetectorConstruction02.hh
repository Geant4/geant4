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
// $Id$
// 
/// \file ExG4DetectorConstruction02.hh
/// \brief Definition of the ExG4DetectorConstruction02 class

#ifndef ExG4DetectorConstruction02_h
#define ExG4DetectorConstruction02_h 1

#include "G4VUserDetectorConstruction.hh"

#include "ExG4DetectorConstruction02Messenger.hh"

#include "G4ThreeVector.hh"

#include "CLHEP/Units/SystemOfUnits.h"

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

/// Simple detector construction with a box volume placed in a world

class ExG4DetectorConstruction02 : public G4VUserDetectorConstruction
{
  public:
    ExG4DetectorConstruction02(
       const G4String& boxMaterialName = "G4_AIR",
       G4double boxHx = 40*CLHEP::cm, 
       G4double boxHy = 40*CLHEP::cm, 
       G4double boxHz = 40*CLHEP::cm,
       const G4String& worldMaterialName = "G4_AIR",
       G4double worldSizeFactor = 1.25);
    ~ExG4DetectorConstruction02();

  public:
    // methods from base class 
    virtual G4VPhysicalVolume* Construct();

    // set methods
    void  SetBoxMaterial(const G4String& materialName);
    void  SetWorldMaterial(const G4String& materialName);
    void  SetBoxDimensions(G4double hx, G4double hy, G4double hz);
    void  SetWorldSizeFactor(G4double factor);
                       
  private:
    ExG4DetectorConstruction02Messenger fMessenger;

    G4String               fBoxMaterialName;
    G4String               fWorldMaterialName;
    G4ThreeVector          fBoxDimensions;
    G4double               fWorldSizeFactor;
    G4LogicalVolume*       fBoxVolume;  
    G4LogicalVolume*       fWorldVolume;  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

