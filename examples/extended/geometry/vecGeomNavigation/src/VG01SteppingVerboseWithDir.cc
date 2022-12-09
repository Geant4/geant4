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
//s
/// \file VG01SteppingVerboseWithDir.cc
/// \brief Implementation of the VG01SteppingVerboseWithDir class

// Author : J. Apostolakis,  (EP/SFT CERN)  2019-21
//
//  Refines the existing G4SteppingVerbose, adding:
//   - direction of particle
//   - safety (at step's start point)

#include "VG01SteppingVerboseWithDir.hh"
#include "G4SteppingManager.hh"

#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

VG01SteppingVerboseWithDir::VG01SteppingVerboseWithDir()
  : G4SteppingVerbose()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

VG01SteppingVerboseWithDir::~VG01SteppingVerboseWithDir()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VG01SteppingVerboseWithDir::Banner()
{
      G4cout << G4endl;
      G4cout << std::setw( 5) << "Step#"     << " "
             << std::setw( 9) << "X (mm)"          << " "
             << std::setw( 9) << "Y (mm)"          << " "
             << std::setw( 9) << "Z (mm)"          << " "
             << std::setw(11) << "Ek (MeV)"     << " "
             << std::setw( 8) << "dE (MeV)"     << " "
             << std::setw( 7) << "Step(mm)"   << " "
             << std::setw(12) << "Track(mm)"   << " "
             << std::setw( 7) << "Safety(mm)" << " | "
             << std::setw(11) << " Dir.x " << " "
             << std::setw(11) << " Dir.y " << " "
             << std::setw(11) << " Dir.z " << " | "
             << std::setw(15) << "Material"
             << std::setw(14) << "Process"
             << std::setw(15) << "NextVolu"
             << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VG01SteppingVerboseWithDir::StepInfo()
{
  CopyState();
  
  G4int prec = G4cout.precision(6);
  
  if( verboseLevel >= 1 ){
    if( verboseLevel >= 4 ) VerboseTrack();
    if( verboseLevel >= 3 ){ Banner(); }

    G4cout << std::setw( 5)<<fTrack->GetCurrentStepNumber() << " "
           << std::setw( 9)<<fTrack->GetPosition().x() / CLHEP::mm << " "
           << std::setw( 9)<<fTrack->GetPosition().y() / CLHEP::mm << " "
           << std::setw( 9)<<fTrack->GetPosition().z() / CLHEP::mm << " " // << " Ek  = "
           << std::setw(11)<<fTrack->GetKineticEnergy() / CLHEP::MeV << " ";
    G4cout.precision(3);
    G4cout << std::setw( 8)<<fStep->GetTotalEnergyDeposit() / CLHEP::MeV << " "
           << std::setw( 7)<<fStep->GetStepLength() / CLHEP::mm << " ";
    G4cout.precision(6);       
    G4cout << std::setw(12)<<fTrack->GetTrackLength() / CLHEP::mm << " ";
    
    G4cout.precision(5);
    G4StepPoint* preStepPt= fTrack->GetStep()->GetPreStepPoint();
    if( preStepPt != nullptr ) {
       G4double safety = preStepPt->GetSafety();
       G4cout  << std::setw( 8) << ( (safety > 1.0e-9) ? safety : 0.0 ) << " ";
    } else {
       G4cout  << std::setw( 8) << "-0.0_(N/A)" << " ";
    }


    G4cout << "| ";
    G4cout << std::setw(11)<<fTrack->GetMomentumDirection().x() << " "
           << std::setw(11)<<fTrack->GetMomentumDirection().y() << " "
           << std::setw(11)<<fTrack->GetMomentumDirection().z() << " | ";
    
    G4cout << std::setw(15) << fTrack->GetMaterial()->GetName() << " ";
    G4VProcess const* pds= fStep->GetPostStepPoint()->GetProcessDefinedStep();
    G4cout << std::setw(14) << ( ( pds != nullptr) ? pds->GetProcessName() :
                                                     G4String( "User Limit") ) <<  " ";
    // if( fStepStatus != fWorldBoundary){     
    if( fTrack->GetNextVolume() != nullptr ) { 
       G4cout << std::setw(10) << fTrack->GetNextVolume()->GetName();
       G4cout << std::setw( 4) << fTrack->GetNextVolume()->GetCopyNo();
    } else {
      G4cout << std::setw(10) << "OutOfWorld";
    }
    
    G4cout << " ";
    
    G4cout << G4endl;
    
    if( verboseLevel == 2 ){
      G4int tN2ndariesTot = fN2ndariesAtRestDoIt +
                            fN2ndariesAlongStepDoIt +
                            fN2ndariesPostStepDoIt;
      if(tN2ndariesTot>0){
        G4cout << "    :----- List of 2ndaries - "
               << "#SpawnInStep=" << std::setw(3) << tN2ndariesTot 
               << "(Rest="  << std::setw(2) << fN2ndariesAtRestDoIt
               << ",Along=" << std::setw(2) << fN2ndariesAlongStepDoIt
               << ",Post="  << std::setw(2) << fN2ndariesPostStepDoIt
               << "), "
               << "#SpawnTotal=" << std::setw(3) << (*fSecondary).size()
               << " ---------------"
               << G4endl;
        
        for(size_t lp1=(*fSecondary).size()-tN2ndariesTot; 
            lp1<(*fSecondary).size(); lp1++){
          G4cout << "    : "
                 << std::setw(8)
                 << (*fSecondary)[lp1]->GetPosition().x() / CLHEP::mm << " "
                 << std::setw(8)
                 << (*fSecondary)[lp1]->GetPosition().y() / CLHEP::mm << " "
                 << std::setw(8)
                 << (*fSecondary)[lp1]->GetPosition().z() / CLHEP::mm << " "
                 << std::setw(8)
                 << (*fSecondary)[lp1]->GetKineticEnergy() / CLHEP::MeV << " "
                 << std::setw(10)
                 << (*fSecondary)[lp1]->GetDefinition()->GetParticleName() << " ";
          G4cout << G4endl;
        }
              
        G4cout << "    :-----------------------------"
               << "----------------------------------"
               << "-- EndOf2ndaries Info ---------------"
               << G4endl;
      }
    }
    
  }
  G4cout.precision(prec);
  // G4cout<< "exit VG01SteppingVerboseWithDir::StepInfo   " <<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VG01SteppingVerboseWithDir::TrackingStarted()
{
  CopyState();
  G4int prec = G4cout.precision(6);
  if( verboseLevel > 0 ){

    Banner();  

    G4cout << std::setw( 5)<<fTrack->GetCurrentStepNumber() << " "
           << std::setw( 9)<<fTrack->GetPosition().x() / CLHEP::mm << " "
           << std::setw( 9)<<fTrack->GetPosition().y() / CLHEP::mm << " "
           << std::setw( 9)<<fTrack->GetPosition().z() / CLHEP::mm << " "  // << " E_0 = "
           << std::setw(11)<<fTrack->GetKineticEnergy() / CLHEP::MeV << " ";
           // << " ( = " << std::setw(8) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy") << " ) "
    G4cout.precision(3);
    G4cout << std::setw( 8) << fStep->GetTotalEnergyDeposit() / CLHEP::MeV << " "
           << std::setw( 7) << fStep->GetStepLength() / CLHEP::mm << " ";
    G4cout.precision(6);
    G4cout << std::setw(12) << fTrack->GetTrackLength() / CLHEP::mm << " ";

    G4Step const* step= fTrack->GetStep();
    G4StepPoint const * preStepPt= nullptr;
    if( step != nullptr )
       preStepPt =fTrack->GetStep()->GetPreStepPoint();
    if( preStepPt != nullptr ) {
       G4double safety = preStepPt->GetSafety();
       G4cout  << std::setw( 8) << ( (safety > 1.0e-9) ? safety : 0.0 ) << " ";
    } else {
       G4cout  << std::setw( 8) << "N/A" << " ";
    }
    G4cout << "| ";
    G4cout << std::setw(11) << fTrack->GetMomentumDirection().x() << " "
           << std::setw(11) << fTrack->GetMomentumDirection().y() << " "
           << std::setw(11) << fTrack->GetMomentumDirection().z() << " | ";

    // Material and volume name
    if(fTrack->GetVolume()){
      G4Material* material = fTrack->GetVolume()->GetLogicalVolume()->GetMaterial();
      if( material != nullptr  ){
         G4cout << std::setw(15) << material->GetName() << " ";
      }
      else {
         G4cout << std::setw(15) << "No-Material" << " ";
      }
      G4cout << std::setw(14) << "initStep" << " ";
      G4cout << std::setw(10) << fTrack->GetVolume()->GetName();
      G4cout << std::setw( 4) << fTrack->GetVolume()->GetCopyNo() << " (initial Volume)";

    } else {
      G4cout << std::setw(10) << "No-Material" << " ";
      G4cout << std::setw(14) << "initStep" << " ";
      G4cout << std::setw(10) << "OutOfWorld" << " ";
    }
    // fTrack->GetMaterial() fails at track start !!!!
    G4cout << G4endl;      
  }
  G4cout.precision(prec);
}
