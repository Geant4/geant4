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
/// \file ExG4DetectorConstruction01.hh
/// \brief Definition of the ExG4DetectorConstruction01 class

#ifndef ExG4DetectorConstruction01_h
#define ExG4DetectorConstruction01_h 1

#include "G4VUserDetectorConstruction.hh"

#include "ExG4DetectorConstruction01Messenger.hh"

#include "G4ThreeVector.hh"

#include "CLHEP/Units/SystemOfUnits.h"

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

/// Simple detector construction with only a world volume

class ExG4DetectorConstruction01 : public G4VUserDetectorConstruction
{
  public:
    ExG4DetectorConstruction01(
       const G4String& materialName = "G4_AIR",
       G4double hx = 50*CLHEP::cm, 
       G4double hy = 50*CLHEP::cm, 
       G4double hz = 50*CLHEP::cm);
    ~ExG4DetectorConstruction01();

  public:
    // methods from base class 
    virtual G4VPhysicalVolume* Construct();

    // set methods
    void  SetMaterial(const G4String& materialName);
    void  SetDimensions(G4double hx, G4double hy, G4double hz);
                       
  private:
    ExG4DetectorConstruction01Messenger fMessenger;

    G4String               fMaterialName;
    G4ThreeVector          fDimensions;
    G4LogicalVolume*       fWorldVolume;  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

