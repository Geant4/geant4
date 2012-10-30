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

#include <iostream>
#include <algorithm>

#include "Tst34EventAction.hh"
#include "Tst34Hit.hh"
#include "G4SystemOfUnits.hh"
#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4UImanager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Event.hh"

#include "G4Timer.hh"
extern G4Timer Timer;
extern G4Timer Timerintern;

Tst34EventAction::Tst34EventAction():
nevent(0),dtime(0.0),calorimeterCollectionId(-1)
{
}

Tst34EventAction::~Tst34EventAction()
{
}

void Tst34EventAction::BeginOfEventAction(const G4Event *evt)
{
  Timerintern.Start();
  G4cout << " ------ Start Tst34EventAction ----- " << G4endl;
  nevent=evt->GetEventID();
  G4cout << " Start generating event Nr " << nevent << G4endl << G4endl;
}

void Tst34EventAction::EndOfEventAction(const G4Event *evt)
{  
  Timerintern.Stop();
  G4cout << G4endl;
  G4cout << "******************************************";
  G4cout << G4endl;
  G4cout << "Internal Real Elapsed Time is: "<< Timerintern.GetRealElapsed();
  G4cout << G4endl;
  G4cout << "Internal System Elapsed Time: " << Timerintern.GetSystemElapsed();
  G4cout << G4endl;
  G4cout << "Internal GetUserElapsed Time: " << Timerintern.GetUserElapsed();
  G4cout << G4endl;
  G4cout << "******************************************" << G4endl;
  dtime+=Timerintern.GetRealElapsed();
  G4cout << " ------ Tst34EventAction::End of event nr. " << nevent
         << "  -----" << G4endl;     

  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4String colNam="Tst34Collection";
  calorimeterCollectionId=SDman->GetCollectionID(colNam);
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  Tst34HitsCollection* THC = 0;
  G4double totE = 0;
  if (calorimeterCollectionId<0) return;

  // Read out of the crysta ECAL

  if (HCE)
  {
    THC=(Tst34HitsCollection *)(HCE->GetHC(calorimeterCollectionId));
  }
  if (THC)
  {
    /// Hits in sensitive Detector
    G4int n_hit = THC->entries();
    G4cout << "  " << n_hit << " hits are stored in Tst34HitsCollection "
           << G4endl;
    G4PrimaryVertex* pvertex=evt->GetPrimaryVertex();   
    ///Computing (x,y,z) of vertex of initial particles  
    G4ThreeVector vtx=pvertex->GetPosition();
    G4PrimaryParticle* pparticle=pvertex->GetPrimary();
    // direction of the Shower
    G4ThreeVector mom=pparticle->GetMomentum()/pparticle->GetMomentum().mag();

    /// For all Hits in sensitive detector
    for (G4int i=0;i<n_hit;i++)
    {
      G4double estep = (*THC)[i]->GetEdep()/GeV;
      if (estep >0.0)
      {
        totE += (*THC)[i]->GetEdep()/GeV;
        G4ThreeVector hitpos=(*THC)[i]->GetPos();
        G4ThreeVector l (hitpos.x(), hitpos.y(), hitpos.z());
        // distance from shower start
      }
    }
  }
  G4cout << " #### Tst34EventAction::Test:"
         << " Total energy deposited in the calorimeter: "
         << totE << " (GeV)" << G4endl;
}
