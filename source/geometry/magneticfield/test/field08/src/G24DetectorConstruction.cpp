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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4FieldManager.hh"
#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4ChordFinder.hh"
#include"G4DormandPrince745.hh"
#include"G4Mag_UsualEqRhs.hh"
#include"G4UserLimits.hh"
#include "G4ClassicalRK4.hh"
#include"G4DormandPrinceRK78.hh"
#include"G4TsitourasRK45.hh"
#include"G4BogackiShampine45.hh"
#include"G4UserLimits.hh"
#include"G4BorisDriver.hh"
#include "G4BorisScheme.hh"
#include"G4BorisDriverSDC.hh"
#include "G4BorisSDC.hh"
class G4VPhysicalVolume;
class G4Material;
class MagneticField;
class G4LogicalVolume;
class G4UserLimits;
#include"G24DetectorConstruction.hpp"

void G24DetectorConstruction::ConstructSDandField()
{
    auto tm = G4TransportationManager::GetTransportationManager();
    G4FieldManager* globalFieldManager = tm->GetFieldManager();
   
    globalFieldManager->SetMinimumEpsilonStep( minEpsilon ) ; 
    globalFieldManager->SetMaximumEpsilonStep( minEpsilon ) ; 
    
    G4MagneticField* magField = new G4UniformMagField(G4ThreeVector(0, 0.0, bz_si*CLHEP::tesla));
    auto fEquation = new G4Mag_UsualEqRhs(magField);
    G4double hmin = 0.01;
    auto stepper = new G4BorisSDC(fEquation);
    auto fdriver = new G4BorisDriverSDC(hmin, stepper);
    auto stepp = new G4BogackiShampine45(fEquation);
    
  auto fChordFinder = new G4ChordFinder( fdriver );
  //auto fChordFinder = new G4ChordFinder( magField,0.01, stepp );
    
    globalFieldManager->SetDetectorField(magField);
    globalFieldManager->SetChordFinder(fChordFinder);
    //globalFieldManager->CreateChordFinder(magField);

}
 G24DetectorConstruction::G24DetectorConstruction()
{}

G24DetectorConstruction::~G24DetectorConstruction()
{}

G4VPhysicalVolume* G24DetectorConstruction::Construct()
{
    G4double world_hx = 100.0*m;
    G4double world_hy = 100.0*m;
    G4double world_hz = 100.0*m;
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");

    G4Box* worldBox
        = new G4Box("World", world_hx, world_hy, world_hz);

    G4LogicalVolume* worldLog
        = new G4LogicalVolume(worldBox, world_mat, "World");
    
    G4VPhysicalVolume* physWorld
        = new G4PVPlacement(0,                     
                            G4ThreeVector(),                   
                            worldLog,              
                            "World",               
                            0,                
                            false,                   
                            0);  

    //setting step length and track length                         
    //G4double circumference = 4*CLHEP::pi/(2.99792458); // for E = 0.2Mev
   
    G4int num_turns =  100;
    G4double max_step = circumference*step_fraction;
    G4double track_len = num_turns*circumference;
    auto fStepLimit = new G4UserLimits(max_step,track_len);
    worldLog->SetUserLimits(fStepLimit);

    return physWorld;                   
}

