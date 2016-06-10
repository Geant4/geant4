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
// $Id: RE05EventAction.cc 66526 2012-12-19 13:41:33Z ihrivnac $
//
/// \file RE05/src/RE05EventAction.cc
/// \brief Implementation of the RE05EventAction class
//

#include "RE05EventAction.hh"

#include "RE05TrackerHit.hh"
#include "RE05CalorimeterHit.hh"
#include "RE05MuonHit.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

RE05EventAction::RE05EventAction()
{
  trackerCollID = -1;
  calorimeterCollID = -1;
  muonCollID = -1;
}

RE05EventAction::~RE05EventAction()
{;}

void RE05EventAction::BeginOfEventAction(const G4Event*)
{
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  if(trackerCollID<0||calorimeterCollID<0||muonCollID<0)
  {
    G4String colNam;
    trackerCollID = SDman->GetCollectionID(colNam="trackerCollection");
    calorimeterCollID = SDman->GetCollectionID(colNam="calCollection");
    muonCollID = SDman->GetCollectionID(colNam="muonCollection");
  }
}

void RE05EventAction::EndOfEventAction(const G4Event* evt)
{
  G4cout << ">>> Event " << evt->GetEventID() << G4endl;
  
  if(trackerCollID<0||calorimeterCollID<0||muonCollID<0) return;

  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  RE05TrackerHitsCollection* THC = 0;
  RE05CalorimeterHitsCollection* CHC = 0;
  RE05MuonHitsCollection* MHC = 0;
  if(HCE)
  {
    THC = (RE05TrackerHitsCollection*)(HCE->GetHC(trackerCollID));
    CHC = (RE05CalorimeterHitsCollection*)(HCE->GetHC(calorimeterCollID));
    MHC = (RE05MuonHitsCollection*)(HCE->GetHC(muonCollID));
  }

  if(THC)
  {
    int n_hit = THC->entries();
    G4cout << "     " << n_hit
         << " hits are stored in RE05TrackerHitsCollection." << G4endl;
  }
  if(CHC)
  {
    int n_hit = CHC->entries();
    G4cout << "     " << n_hit
         << " hits are stored in RE05CalorimeterHitsCollection." << G4endl;
    G4double totE = 0;
    for(int i=0;i<n_hit;i++)
    { totE += (*CHC)[i]->GetEdep(); }
    G4cout << "     Total energy deposition in calorimeter : "
         << totE / GeV << " (GeV)" << G4endl;
  }
  if(MHC)
  {
    int n_hit = MHC->entries();
    G4cout << "     " << n_hit
         << " hits are stored in RE05MuonHitsCollection." << G4endl;
  }
}

