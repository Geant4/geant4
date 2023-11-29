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
/// \file XraySPODetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// Authors: P.Dondero (paolo.dondero@cern.ch), R.Stanzani (ronny.stanzani@cern.ch)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#ifndef XraySPODetectorConstruction_h
#define XraySPODetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Cache.hh"
#include "G4GDMLParser.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class XraySPODetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class XraySPODetectorConstruction : public G4VUserDetectorConstruction
{
    public:

      explicit XraySPODetectorConstruction();
      ~XraySPODetectorConstruction() override = default;
      G4VPhysicalVolume* Construct() override;
      G4VPhysicalVolume* ConstructDetector();
      void SetReadFile(G4String &);
      G4VPhysicalVolume* ConstructVolumes();
      void DefineMaterials();

    private:

      G4GDMLParser fParser;
      G4String fReadFile;
      G4VPhysicalVolume* fPhysiWorld;
      G4Material* fDefaultMaterial;
      XraySPODetectorMessenger* fDetectorMessenger;
      G4LogicalVolume* fWorld_log;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
