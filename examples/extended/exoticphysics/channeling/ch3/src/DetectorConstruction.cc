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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4RegionStore.hh"
#include "G4VisAttributes.hh"
#include "PrimaryGeneratorAction.hh"
#include <CLHEP/Units/SystemOfUnits.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    //Check overlap option
    G4bool checkOverlaps = true;

    //Materials
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
    G4Material* silicon = nist->FindOrBuildMaterial("G4_Si");

    //World
    G4Box* solidWorld = new G4Box("World", 0.2*CLHEP::m, 0.2*CLHEP::m, 0.3*CLHEP::m);
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");  
    G4VPhysicalVolume* physWorld = new G4PVPlacement
                                    (0,                    // no rotation
                                    G4ThreeVector(),       // centre position
                                    logicWorld,            // its logical volume
                                    "World",               // its name
                                    0,                     // its mother volume
                                    false,                 // no boolean operation
                                    0,                     // copy number
                                    checkOverlaps);        // overlaps checking
    logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());


    // --------------- Crystal ------------------------------------

    //Select crystal material
    fCrystalMaterial = nist->FindOrBuildMaterial("G4_W");

    //Setting crystal rotation angle (also the angle of crystal planes vs the beam)
    //Crystal rotation angle (also the angle of crystal planes vs the beam)
    G4double angleX = 0.*1e-6; //rad
    G4RotationMatrix* crystalRotationMatrix = new G4RotationMatrix;
    crystalRotationMatrix->rotateY(-angleX);

    //setting crystal dimensions
    /*at high energies the electromagnetic shower in
      oriented tungsten should behave similarly to the
      e.m. shower in several times thicker amorphous tungsten
      (several times reduction of the effective radiation length)*/
    G4ThreeVector crystalSize = G4ThreeVector(10.*CLHEP::mm,
                                              10.*CLHEP::mm,
                                              0.1*CLHEP::mm);

    //Setting crystal position
    G4ThreeVector posCrystal = G4ThreeVector(0., 0., crystalSize.z()/2.);

    //crystal volume
    G4Box* solidCrystal = new G4Box("Crystal",
                                    crystalSize.x()/2,
                                    crystalSize.y()/2,
                                    crystalSize.z()/2.);
    
    fLogicCrystal = new G4LogicalVolume(solidCrystal,
                                       fCrystalMaterial,
                                       "Crystal");
    new G4PVPlacement(crystalRotationMatrix,
                      posCrystal,
                      fLogicCrystal,
                      "Crystal",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    
    //crystal region (necessary for the FastSim model)
    G4Region* regionCh = new G4Region("Crystal");
    regionCh->AddRootLogicalVolume(fLogicCrystal);

    //visualization attributes
    G4VisAttributes* crystalVisAttribute =
        new G4VisAttributes(G4Colour(1., 0., 0.));
    crystalVisAttribute->SetForceSolid(true);
    fLogicCrystal->SetVisAttributes(crystalVisAttribute);
        
    //print crystal info
    G4cout << "Crystal size: " << crystalSize.x()/CLHEP::mm
           << " " << crystalSize.y()/CLHEP::mm
           << " " << crystalSize.z()/CLHEP::mm << " mm3" << G4endl;
    G4cout << "Crystal angleX: " << angleX << " rad" << G4endl;

    // --------------- Detector -----------------------------------
    //Setting detector position
    G4ThreeVector posDetector = G4ThreeVector(0, 0, 0.1*CLHEP::m);

    //particle detector volume
    G4Box* detector = new G4Box("Detector",
                                10*CLHEP::cm/2,
                                10*CLHEP::cm/2,
                                0.3*CLHEP::mm/2);
    
    G4LogicalVolume* logicDetector = new G4LogicalVolume(detector,
                                                         silicon,
                                                         "Detector");
    new G4PVPlacement(0, 
                      posDetector, 
                      logicDetector,
                      "Detector", 
                      logicWorld, 
                      false, 
                      0, 
                      checkOverlaps);

    //visualization attributes
    G4VisAttributes* detectorVisAttribute =
        new G4VisAttributes(G4Colour(0., 0., 1));
    detectorVisAttribute->SetForceSolid(true);
    logicDetector->SetVisAttributes(detectorVisAttribute);
    
    //always return the physical World
    return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
    // --------------- fast simulation ----------------------------
    //extract the region of the crystal from the store
    G4RegionStore* regionStore = G4RegionStore::GetInstance();
    G4Region* regionCh = regionStore->GetRegion("Crystal");

    //create the channeling model for this region
    G4ChannelingFastSimModel* channelingModel =
        new G4ChannelingFastSimModel("ChannelingModel", regionCh);

    //Crystal planes or axes considered
    //Use brackets (...) for planes and <...> for axes
    G4String lattice = "<111>";

    //activate the channeling model
    channelingModel->Input(fCrystalMaterial, lattice);

    /*
    activate radiation model (do it only when you want to take into account the
    radiation production in an oriented crystal; it reduces simulation speed.)
    */
    G4bool activateRadiationModel = true;
    if (activateRadiationModel)
    {
        channelingModel->RadiationModelActivate();
        G4cout << "Radiation model activated" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

