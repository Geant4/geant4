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
/// \file Par02DetectorConstruction.cc
/// \brief Implementation of the Par02DetectorConstruction class

#include "Par02DetectorConstruction.hh"
#include "G4ProductionCuts.hh"
#include "G4SystemOfUnits.hh"
#include "G4RegionStore.hh"
#include "G4GDMLParser.hh"
#include "G4AutoDelete.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02DetectorConstruction::Par02DetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02DetectorConstruction::~Par02DetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* Par02DetectorConstruction::Construct() {
  G4GDMLParser parser;
  parser.Read( "Par02FullDetector.gdml" );
  G4cout << "Geometry loaded from  file .......Par02FullDetector.gdml " << G4endl;

  // This GDML detector description uses the auxiliary information part to store
  // information regarding which Geant4 volumes have a fast simulation model.

  const G4GDMLAuxMapType* aAuxMap = parser.GetAuxMap();
  for ( G4GDMLAuxMapType::const_iterator iter = aAuxMap->begin();
        iter != aAuxMap->end(); ++iter ) {
    for ( G4GDMLAuxListType::const_iterator vit = (*iter).second.begin();
          vit != (*iter).second.end(); ++vit ) {
      if ( (*vit).type == "FastSimModel" ) {
        G4LogicalVolume* myvol = (*iter).first;
        if ( ( myvol->GetName() ).find( "Tracker" ) != std::string::npos ) {
          fTrackerList.push_back( new G4Region( myvol->GetName() ) );
          fTrackerList.back()->AddRootLogicalVolume( myvol );
          G4cout << G4endl << "tracker !!!" << G4endl;
        } else if ( ( myvol->GetName() ).find( "HCal" ) != std::string::npos ) {
          fHCalList.push_back( new G4Region( myvol->GetName() ) );
          fHCalList.back()->AddRootLogicalVolume( myvol );
          G4cout << G4endl << "hcal !!!" << G4endl;
        } else if ( ( myvol->GetName() ).find( "ECal" ) != std::string::npos ) {
          fECalList.push_back( new G4Region( myvol->GetName() ) );
          fECalList.back()->AddRootLogicalVolume( myvol );
          G4cout << G4endl << "ecal !!!" << G4endl;
        } else if ( ( myvol->GetName() ).find( "Muon" ) != std::string::npos ) {
          fMuonList.push_back( new G4Region( myvol->GetName() ) );
          fMuonList.back()->AddRootLogicalVolume( myvol );
        } else {
          G4cout << G4endl << "NOT A KNOWN DETECTOR !!!" << G4endl;
        }
      }
    }
  }
  for ( G4int iterTracker = 0; iterTracker < G4int( fTrackerList.size() ); 
        iterTracker++ ) {
    fTrackerList[ iterTracker ]->SetProductionCuts( new G4ProductionCuts() );
    fTrackerList[ iterTracker ]->GetProductionCuts()->SetProductionCut
      ( 1.0* ( ( *fTrackerList[ iterTracker ]->GetRootLogicalVolumeIterator() )->
          GetMaterial()->GetRadlen() ) );
    fTrackerList[ iterTracker ]->GetProductionCuts()->
      SetProductionCut( 1.0*m, idxG4GammaCut );
  }
  for ( G4int iterECal = 0; iterECal < G4int( fECalList.size() ); iterECal++ ) {
    fECalList[ iterECal ]->SetProductionCuts( new G4ProductionCuts() );
    fECalList[ iterECal ]->GetProductionCuts()->SetProductionCut
      ( 0.5* ( ( *fECalList[ iterECal ]->GetRootLogicalVolumeIterator() )->
          GetMaterial()->GetRadlen() ) );
    fECalList[ iterECal ]->GetProductionCuts()->
      SetProductionCut( 0.1*m, idxG4GammaCut );
  }
  for ( G4int iterHCal = 0; iterHCal < G4int( fHCalList.size() ); iterHCal++ ) {
    fHCalList[ iterHCal ]->SetProductionCuts( new G4ProductionCuts() );
    fHCalList[ iterHCal ]->GetProductionCuts()->SetProductionCut(
      0.5* ( ( *fHCalList[iterHCal]->GetRootLogicalVolumeIterator() )->
        GetMaterial()->GetRadlen() ) );
    fHCalList[ iterHCal ]->GetProductionCuts()->
      SetProductionCut( 1.0*m, idxG4GammaCut );
  }

  // Returns the pointer to the physical world.
  return parser.GetWorldVolume();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par02DetectorConstruction::ConstructSDandField() {
  for ( G4int iterTracker = 0; iterTracker < G4int( fTrackerList.size() ); 
        iterTracker++ ) {
    // Bound the fast simulation model for the tracker subdetector
    // to all the corresponding Geant4 regions
    Par02FastSimModelTracker* fastSimModelTracker
      = new Par02FastSimModelTracker( "fastSimModelTracker", fTrackerList[ iterTracker ],
                                      Par02DetectorParametrisation::eCMS );
                                     
    // Register the fast simulation model for deleting
    G4AutoDelete::Register(fastSimModelTracker);
  }
  for ( G4int iterECal = 0; iterECal < G4int( fECalList.size() ); iterECal++ ) {
    // Bound the fast simulation model for the electromagnetic calorimeter
    // to all the corresponding Geant4 regions
    Par02FastSimModelEMCal* fastSimModelEMCal
      = new Par02FastSimModelEMCal( "fastSimModelEMCal", fECalList[ iterECal ],
                                    Par02DetectorParametrisation::eCMS );
                                     
    // Register the fast simulation model for deleting
    G4AutoDelete::Register(fastSimModelEMCal);
  }
  for ( G4int iterHCal = 0; iterHCal < G4int( fHCalList.size() ); iterHCal++ ) {
    // Bound the fast simulation model for the hadronic calorimeter
    // to all the corresponding Geant4 regions
    Par02FastSimModelHCal* fastSimModelHCal
      = new Par02FastSimModelHCal( "fastSimModelHCal", fHCalList[ iterHCal ],
                                   Par02DetectorParametrisation::eCMS );
                                     
    // Register the fast simulation model for deleting
    G4AutoDelete::Register( fastSimModelHCal );
  }
  // Currently we don't have a fast muon simulation model to be bound
  // to all the corresponding Geant4 regions.
  // But it could be added in future, in a similar way as done above for
  // the tracker subdetector and the electromagnetic and hadronic calorimeters.

  // Add global magnetic field
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger( fieldValue );
  fMagFieldMessenger->SetVerboseLevel(1);
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

