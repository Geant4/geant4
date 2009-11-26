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
//

#include "WLSEventAction.hh"

#include "WLSRunAction.hh"

#include "WLSEventActionMessenger.hh"

#include "WLSPhotonDetHit.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "WLSTrajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"

//#include "G4ThreeVector.hh"

//#include "G4UnitsTable.hh"

#include "Randomize.hh"

// Purpose: Invoke visualization at the end
//          Also can accumulate statistics regarding hits
//          in the PhotonDet detector

WLSEventAction::WLSEventAction(WLSRunAction* RA)
 : runaction(RA), verboselevel(0),
   printModulo(100), drawFlag("all")
{
  eventMessenger = new WLSEventActionMessenger(this);

  forcedrawphotons = false;
  forcenophotons   = false;

}

WLSEventAction::~WLSEventAction()
{
  delete eventMessenger;
}

void WLSEventAction::BeginOfEventAction(const G4Event* evt)
{
 G4int evtNb = evt->GetEventID();
 if (evtNb%printModulo == 0)
    G4cout << "\n---> Begin of Event: " << evtNb << G4endl;

 if(verboselevel>0)
    G4cout << "<<< Event  " << evtNb << " started." << G4endl;
}

void WLSEventAction::EndOfEventAction(const G4Event* evt)
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  // Visualization of Trajectory
  if(pVVisManager)
  {
   G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();

   G4int n_trajectories = 0;
   if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
   G4cout << "n_trajectories: " << n_trajectories << G4endl;
   if (drawFlag == "all") G4cout << "draw all trajectories" << G4endl;
   if (drawFlag == "charged") G4cout << "draw only charged" << G4endl;

   for(G4int i=0; i<n_trajectories; i++)
      { WLSTrajectory* trj = 
                      (WLSTrajectory *)((*(evt->GetTrajectoryContainer()))[i]);
        if (drawFlag == "all") {
           G4cout << "Now calling DrawTrajectory" << G4endl;
           G4cout << "Particle Name: " << trj->GetParticleName() << G4endl;
           trj->DrawTrajectory(50);
        }
        else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
                               trj->DrawTrajectory(50);
        else if (trj->GetParticleName()=="opticalphoton")
        {
          G4cout << "We should be drawing an opticalphoton" << G4endl;
          trj->SetForceDrawTrajectory(forcedrawphotons);
          trj->SetForceNoDrawTrajectory(forcenophotons);
          trj->DrawTrajectory(50);
        }
      }
  }

  if (verboselevel>0)
     G4cout << "<<< Event  " << evt->GetEventID() << " ended." << G4endl;
 
  // Save the Random Engine Status at the of event
  if (runaction->GetRndmFreq() == 2)
  {
     CLHEP::HepRandom::saveEngineStatus("endOfEvent.rndm");
     G4int evtNb = evt->GetEventID();
     if (evtNb%printModulo == 0)
     {
        G4cout << "\n---> End of Event: " << evtNb << G4endl;
        CLHEP::HepRandom::showEngineStatus();
     }
  }

  // Get Hits from the detector if any
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4String colName = "PhotonDetHitCollection";
  mppcCollID = SDman->GetCollectionID(colName);

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  WLSPhotonDetHitsCollection* mppcHC = 0;

  // Get the hit collections
  if (HCE)
  {
     if (mppcCollID>=0) mppcHC = 
                        (WLSPhotonDetHitsCollection*)(HCE->GetHC(mppcCollID));
  }

  // Get hit information about photons that reached the detector in this event
  if (mppcHC)
  {
//     G4int n_hit = mppcHC->entries();
  }
}

G4int WLSEventAction::GetEventNo()
{
  G4int evno = fpEventManager->GetConstCurrentEvent()->GetEventID();
  return evno;
}

void WLSEventAction::SetEventVerbose(G4int level)
{
  verboselevel = level;
}
