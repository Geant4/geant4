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
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// If you use this example, please cite the following publication:
// Rad. Prot. Dos. 133 (2009) 2-11

#include "G4Event.hh"
#include "G4AnalysisManager.hh"
#include "Randomize.hh"

#include "EventAction.hh"
#include "RunAction.hh"

EventAction::EventAction(RunAction* run)
:fRun(run)
{}

EventAction::~EventAction()
{}

void EventAction::BeginOfEventAction(const G4Event* evt)
{  
  G4int evtNb = evt->GetEventID();
  fRun->SetNumEvent(evtNb);
  fRun->SetDoseN(0);
  fRun->SetDoseC(0);
}

void EventAction::EndOfEventAction(const G4Event* )
{  
  G4AnalysisManager* man = G4AnalysisManager::Instance();

  // Save total absorbed dose in phantom

  if (fRun->GetDoseN()>0 || fRun->GetDoseC()>0) 
  {
    // Fill ntuple #4
    man->FillNtupleDColumn(4,0,fRun->GetDoseN());
    man->FillNtupleDColumn(4,1,fRun->GetDoseC());
    man->AddNtupleRow(4);

    G4cout << "   ===> The incident alpha particle has reached the targeted cell :" << G4endl;
    G4cout << "   -----> total absorbed dose within Nucleus   is (Gy) = " << fRun->GetDoseN() << G4endl;
    G4cout << "   -----> total absorbed dose within Cytoplasm is (Gy) = " << fRun->GetDoseC() << G4endl;
    G4cout << G4endl;
  }
  else
  {
    G4cout << "   ===> Sorry, the incident alpha particle has missed the targeted cell !" << G4endl;
    G4cout << G4endl;
  }
}
