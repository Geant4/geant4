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
// $Id: RE01EventAction.cc,v 1.2 2006-06-29 17:43:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//


#include "RE01EventAction.hh"

#include "RE01TrackerHit.hh"
#include "RE01CalorimeterHit.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "RE01Trajectory.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"

RE01EventAction::RE01EventAction()
{
  trackerCollID = -1;
  calorimeterCollID = -1;
  muonCollID = -1;
}

RE01EventAction::~RE01EventAction()
{;}

void RE01EventAction::BeginOfEventAction(const G4Event*)
{
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  if(trackerCollID<0||calorimeterCollID<0||muonCollID<0)
  {
    G4String colNam;
    trackerCollID = SDman->GetCollectionID(colNam="trackerCollection");
    calorimeterCollID = SDman->GetCollectionID(colNam="calCollection");
  }
}

void RE01EventAction::EndOfEventAction(const G4Event* evt)
{
  G4cout << ">>> Summary of Event " << evt->GetEventID() << G4endl;
  
  if(trackerCollID<0||calorimeterCollID<0) return;

  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  RE01TrackerHitsCollection* THC = 0;
  RE01CalorimeterHitsCollection* CHC = 0;
  if(HCE)
  {
    THC = (RE01TrackerHitsCollection*)(HCE->GetHC(trackerCollID));
    CHC = (RE01CalorimeterHitsCollection*)(HCE->GetHC(calorimeterCollID));
  }

  if(THC)
  {
    int n_hit = THC->entries();
    G4cout << G4endl;
    G4cout << "Tracker hits --------------------------------------------------------------" << G4endl;
    G4cout << n_hit << " hits are stored in RE01TrackerHitsCollection." << G4endl;
    G4cout << "List of hits in tracker" << G4endl;
    for(int i=0;i<n_hit;i++)
    { (*THC)[i]->Print(); }
  }
  if(CHC)
  {
    int n_hit = CHC->entries();
    G4cout << G4endl;
    G4cout << "Calorimeter hits --------------------------------------------------------------" << G4endl;
    G4cout << n_hit << " hits are stored in RE01CalorimeterHitsCollection." << G4endl;
    G4double totE = 0;
    for(int i=0;i<n_hit;i++)
    { totE += (*CHC)[i]->GetEdep(); }
    G4cout << "     Total energy deposition in calorimeter : "
         << totE / GeV << " (GeV)" << G4endl;
  }

  // get number of stored trajectories
  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  // extract the trajectories and print them out
  G4cout << G4endl;
  G4cout << "Trajectories in tracker --------------------------------------------------------------" << G4endl;
  for(G4int i=0; i<n_trajectories; i++) 
  {
    RE01Trajectory* trj = (RE01Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
    trj->ShowTrajectory();
  }
    
  G4cout << G4endl;
  G4cout << "Primary particles --------------------------------------------------------------" << G4endl;
  G4int n_vertex = evt->GetNumberOfPrimaryVertex();
  for(G4int iv=0;iv<n_vertex;iv++)
  {
    G4PrimaryVertex* pv = evt->GetPrimaryVertex(iv);
    G4cout << G4endl;
    G4cout << "Primary vertex "
           << G4ThreeVector(pv->GetX0(),pv->GetY0(),pv->GetZ0())
           << "   at t = " << (pv->GetT0())/ns << " [ns]" << G4endl;
    G4PrimaryParticle* pp = pv->GetPrimary();
    while(pp)
    {
      PrintPrimary(pp,0);
      pp = pp->GetNext();
    }
  }
}

void RE01EventAction::PrintPrimary(G4PrimaryParticle* pp,G4int ind)
{
  for(G4int ii=0;ii<=ind;ii++)
  { G4cout << "  "; }
  G4cout << "==PDGcode " << pp->GetPDGcode() << " ";
  if(pp->GetG4code()!=0)
  { G4cout << "(" << pp->GetG4code()->GetParticleName() << ")"; }
  else
  { G4cout << "is not defined in G4"; }
  G4cout << " " << pp->GetMomentum()/GeV << " [GeV] ";
  if(pp->GetTrackID()<0)
  { G4cout << G4endl; }
  else
  { G4cout << ">>> G4Track ID " << pp->GetTrackID() << G4endl; }

  G4PrimaryParticle* daughter = pp->GetDaughter();
  while(daughter)
  {
    PrintPrimary(daughter,ind+1);
    daughter = daughter->GetNext();
  }
}


