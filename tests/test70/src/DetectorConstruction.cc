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
#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"

#ifdef G4VIS_USE
#include "G4VisAttributes.hh"
#endif

//______________________________________________________________________
DetectorConstruction::DetectorConstruction()
{
    PWorld = NULL;
    LWorld = NULL;
    SWorld = NULL;
    waterMaterial = NULL;
}

//______________________________________________________________________
DetectorConstruction::~DetectorConstruction()
{
    if(PWorld)
        delete PWorld ;
}

//______________________________________________________________________
G4VPhysicalVolume* DetectorConstruction::Construct()

{
    DefineMaterials();
    return ConstructDetector();
}

//______________________________________________________________________
void DetectorConstruction::DefineMaterials()
{

    // Water is defined from NIST material database
    G4NistManager * man = G4NistManager::Instance();
    G4Material * H2O = man->FindOrBuildMaterial("G4_WATER");

    // Default materials in setup.
    waterMaterial = H2O;

}

//______________________________________________________________________
G4VPhysicalVolume* DetectorConstruction::ConstructDetector()
{
    double  WorldSizeX  = 10.*cm ;
    double  WorldSizeY  = 10.*cm ;
    double  WorldSizeZ  = 10.*cm ;
    SWorld = new G4Box("DefaultRegionForTheWorld", WorldSizeX,WorldSizeY,WorldSizeZ);
    LWorld = new G4LogicalVolume(SWorld,waterMaterial,"DefaultRegionForTheWorld");

    PWorld = new G4PVPlacement(0,                           //no rotation
                               G4ThreeVector(),             //at (0,0,0)
                               "DefaultRegionForTheWorld",  //its name
                               LWorld,                      //its logical volume
                               NULL,                        //its mother  volume
                               false,                       //no boolean operation
                               0);                          //copy number

    //__________________________________________________________________
    // Visualization attributes
#ifdef G4VIS_USE
    G4VisAttributes* worldVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0)); //White
    worldVisAtt->SetVisibility(true);
    LWorld->SetVisAttributes(worldVisAtt);
#endif

    return PWorld;
}
