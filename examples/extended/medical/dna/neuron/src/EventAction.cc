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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// and papers
// M. Batmunkh et al. J Radiat Res Appl Sci 8 (2015) 498-507
// O. Belov et al. Physica Medica 32 (2016) 1510-1520
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// -------------------------------------------------------------------
// November 2016
// -------------------------------------------------------------------
//
// $Id$
//
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "G4Event.hh"
#include "Randomize.hh"
#include "EventAction.hh"
#include "RunAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include "Analysis.hh"
#include "Run.hh"
#include "G4RunManager.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::EventAction(RunAction* run)
:fRunAction(run)
{
 remove ("OutputPerEvent.out");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::~EventAction()
{
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::BeginOfEventAction(const G4Event* evt)
{  
  G4int evtNb = evt->GetEventID();
  fRunAction->SetNumEvent(evtNb);
  // 
  fRunAction->SetEdepALL(0);
  fRunAction->SetEdepMedium(0);
  fRunAction->SetEdepSlice(0);
  fRunAction->SetEdepNeuron(0);
  fRunAction->SetEdepSoma(0);
  fRunAction->SetEdepDend(0);
  fRunAction->SetEdepAxon(0);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::EndOfEventAction(const G4Event* evt)
{  
 //Get Postion and Momentum of primary
 // beam index of per track
 G4int evtNb = evt->GetEventID();
 // const G4ThreeVector& pos = evt->GetPrimaryVertex()->GetPosition();
 // const G4ThreeVector& mom = evt->GetPrimaryVertex()
 //                               ->GetPrimary()->GetMomentum();

 Run* run = static_cast<Run*>(
            G4RunManager::GetRunManager()->GetNonConstCurrentRun());
 run->AddEdepALL(fRunAction->GetEdepALL()); 
 run->AddEdepMedium(fRunAction->GetEdepMedium()); 
 run->AddEdepSlice(fRunAction->GetEdepSlice()); 
 run->AddEdepNeuron(fRunAction->GetEdepNeuron()); 
 run->AddEdepSoma(fRunAction->GetEdepSoma()); 
 run->AddEdepDend(fRunAction->GetEdepDend()); 
 run->AddEdepAxon(fRunAction->GetEdepAxon()); 
   // to calculate LET
   /*
   if (fRunAction->GetEdepALL() > 0.) 
   {   
   G4double LET = fRunAction->GetEdepALL()
       / ftrackLength() ;  
   run->AddPrimaryLET(LET);
   }   
   */
/*    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleDColumn(4,0,fRunAction->GetEdepALL()/keV );
    analysisManager->FillNtupleDColumn(4,1,fRunAction->GetEdepMedium()/keV );
    analysisManager->FillNtupleDColumn(4,2,
           (fRunAction->GetEdepSlice()+fRunAction->GetEdepNeuron())/keV);
    analysisManager->FillNtupleDColumn(4,3,fRunAction->GetEdepNeuron()/keV );
    analysisManager->AddNtupleRow(4);
*/
  std::ofstream OutputEdep("OutputPerEvent.out", std::ios::app);
  OutputEdep<<   evtNb+1          << '\t' << "   "  // event number
        // edep in all volume
   <<   fRunAction->GetEdepALL()/keV         << '\t' << "   " 
        // outside bounding box 
   <<   fRunAction->GetEdepMedium()/keV         << '\t' << "   " 
   <<   (fRunAction->GetEdepSlice()+fRunAction->GetEdepNeuron())/keV 
   << '\t' << "   "//  inside Bounding Slice
   <<   fRunAction->GetEdepNeuron()/keV         << '\t' << "   " 
   <<   fRunAction->GetEdepSoma()/keV         << '\t' << "   " 
   <<   fRunAction->GetEdepDend()/keV         << '\t' << "   " 
   <<   fRunAction->GetEdepAxon()/keV         << '\t' << "   " 
   << G4endl;  
 
 /*
 if (fRunAction->GetEdepNeuron() > 0.)
 {
 std::ofstream WriteDataInNeurN("0-Neur-EventN.out", std::ios::app);
 WriteDataInNeurN <<   evtNb+1          << '\t' << "   " // 
   <<   pos.x()/um          << '\t' << "   "  
   <<   pos.y()/um          << '\t' << "   "  
   <<   pos.z()/um          << '\t' << "   " 
   << G4endl; 
 }  
   */
}
