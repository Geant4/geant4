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
/// \file field/field05/src/F05SteppingVerbose.cc
/// \brief Implementation of the F05SteppingVerbose class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#include "F05SteppingVerbose.hh"

#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

F05SteppingVerbose::F05SteppingVerbose()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

F05SteppingVerbose::~F05SteppingVerbose()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void F05SteppingVerbose::StepInfo()
{
  CopyState();
  
  G4int prec = G4cout.precision(6);

  if( verboseLevel >= 1 ){
    if( verboseLevel >= 4 ) VerboseTrack();
    if( verboseLevel >= 3 ){
      G4cout << G4endl;    
      G4cout << std::setw( 5) << "#Step#"     << " "
             << std::setw(10) << "X"          << "    "
             << std::setw(10) << "Y"          << "    "  
             << std::setw(10) << "Z"          << "    "
             << std::setw(10) << "Direction x"      << " "
             << std::setw(10) << "dir y"      << " "  
             << std::setw(10) << "dir z"      << " "
             << std::setw( 8) << "Time"      << " "
             << std::setw( 6) << "dEStep"     << " "  
             << std::setw(10) << "StepLeng"     
             << std::setw(10) << "TrakLeng" 
             << std::setw(10) << "NextVolu" 
             << std::setw(10) << "Process"    << " "
             << G4endl;                  
    }

    G4cout << std::setw( 5) << fTrack->GetCurrentStepNumber() << " "
           << std::setw(10) << G4BestUnit(fTrack->GetPosition().x(),"Length")
           << std::setw(10) << G4BestUnit(fTrack->GetPosition().y(),"Length")
           << std::setw(10) << G4BestUnit(fTrack->GetPosition().z(),"Length");
    G4cout << std::setw(10) << fTrack->GetMomentumDirection().x() << " "
           << std::setw(10) << fTrack->GetMomentumDirection().y() << " "
           << std::setw(10) << fTrack->GetMomentumDirection().z() << " ";
    G4cout << std::setw( 8) << G4BestUnit(fTrack->GetGlobalTime(),"Time")
           << std::setw( 6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
           << std::setw(10) << G4BestUnit(fStep->GetStepLength(),"Length")
           << std::setw(10) << G4BestUnit(fTrack->GetTrackLength(),"Length");

    if( fTrack->GetNextVolume() != 0 ) { 
      G4cout << std::setw(10) << fTrack->GetVolume()->GetName();
    } else {
      G4cout << std::setw(10) << "OutOfWorld";
    }

    if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != 0){
      G4cout << "  " 
             << std::setw(10) << fStep->GetPostStepPoint()->GetProcessDefinedStep()
                                        ->GetProcessName();
    } else {
      G4cout << "   UserLimit";
    }
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
                   lp1<(*fSecondary).size(); lp1++)
        {
          G4cout << "    : "
                 << std::setw(6)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
                 << std::setw(6)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
                 << std::setw(6)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
                 << std::setw(6)
                 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
                 << std::setw(10)
                 << (*fSecondary)[lp1]->GetDefinition()->GetParticleName();
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void F05SteppingVerbose::TrackingStarted()
{

  CopyState();
  G4int prec = G4cout.precision(6);
  if( verboseLevel > 0 ){

    G4cout << std::setw( 5) << "Step#"      << " "
           << std::setw(10) << "X"          << "    "
           << std::setw(10) << "Y"          << "    "  
           << std::setw(10) << "Z"          << "    "
           << std::setw(10) << "Direction x"  << " "
           << std::setw(10) << "dir y"        << " "  
           << std::setw(10) << "dir z"        << " "
           << std::setw(10) << "Time"      << " "
           << std::setw(10) << "dEStep"     << " "  
           << std::setw(10) << "StepLeng"  
           << std::setw(10) << "TrakLeng"
           << std::setw(10) << "Volume"     << "  "
           << std::setw(10) << "Process"    
           << G4endl;             

    G4cout << std::setw( 5) << fTrack->GetCurrentStepNumber() << " "
           << std::setw(10) << G4BestUnit(fTrack->GetPosition().x(),"Length")
           << std::setw(10) << G4BestUnit(fTrack->GetPosition().y(),"Length")
           << std::setw(10) << G4BestUnit(fTrack->GetPosition().z(),"Length");
    G4cout << std::setw(10) << fTrack->GetMomentumDirection().x() << " "
           << std::setw(10) << fTrack->GetMomentumDirection().y() << " "
           << std::setw(10) << fTrack->GetMomentumDirection().z() << " ";
    G4cout << std::setw( 8) << G4BestUnit(fTrack->GetGlobalTime(),"Time")
           << std::setw( 6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
           << std::setw(10) << G4BestUnit(fStep->GetStepLength(),"Length")
           << std::setw(10) << G4BestUnit(fTrack->GetTrackLength(),"Length");

    if(fTrack->GetNextVolume()){
      G4cout << std::setw(10) << fTrack->GetVolume()->GetName();
    } else {
      G4cout << std::setw(10) << "OutOfWorld";
    }
    G4cout  << "    initStep" << G4endl;
  }
  G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
