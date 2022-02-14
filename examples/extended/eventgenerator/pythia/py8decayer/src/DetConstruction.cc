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
/// \file eventgenerator/pythia/pythia8decayer/src/DetConstruction.cc
/// \brief Implementation of the DetConstruction class
///
/// \author J. Yarba; FNAL

#include "DetConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetConstruction::Construct()
{

   // Define materials via NIST manager
  
   auto worldMaterial = 
      G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"); 
   auto detMaterial   = 
      G4NistManager::Instance()->FindOrBuildMaterial("G4_C"); 

   // World volume
 
   G4ThreeVector worldSize( 200.*CLHEP::cm, 200.*CLHEP::cm, 200.*CLHEP::cm );
  

   fWorld = new G4PVPlacement( 
      0, G4ThreeVector(), 
      "world_phys",
      new G4LogicalVolume( new G4Box( "world_solid", 
                                      0.5*worldSize.x(), 
                                      0.5*worldSize.y(), 
                                      0.5*worldSize.z() ),
                           worldMaterial, "world_log", 0, 0, 0 ), 
      0,
      false, 0 );

  // "Detector"

  double rmin    = 10.*CLHEP::cm;
  double rmax    = 80.*CLHEP::cm;
  double zlength = 98.*CLHEP::cm;
  
  G4VSolid* sDet = new G4Tubs( "det_solid", 
                               rmin, rmax, 0.5*zlength, 
                               0., 2.*CLHEP::pi );
  
  G4LogicalVolume* lDet = new G4LogicalVolume( sDet, detMaterial, "det_log", 
                                               0, 0, 0 );
  
  fDet = new G4PVPlacement( 0, G4ThreeVector(), "det_phys",
                            lDet, 
                            fWorld,  // it's mother (physical) volume 
                            false, 0 );
  
  
  //always return the root volume
  //
  return fWorld;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

