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
/// \file field/field04/src/F04EventAction.cc
/// \brief Implementation of the F04EventAction class
//
//
#include "F04EventAction.hh"

#include "F04RunAction.hh"

#include "F04EventActionMessenger.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VVisManager.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04EventAction::F04EventAction(F04RunAction* runAction)
 : fRunaction(runAction), fVerboselevel(0),
   fPrintModulo(10), fDrawFlag("all")
{
  fEventMessenger = new F04EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04EventAction::~F04EventAction()
{
  delete fEventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04EventAction::BeginOfEventAction(const G4Event* evt)
{
 G4int evtNb = evt->GetEventID();
 if (evtNb%fPrintModulo == 0)
    G4cout << "\n---> Begin of Event: " << evtNb << G4endl;

 if(fVerboselevel>0)
    G4cout << "<<< Event  " << evtNb << " started." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04EventAction::EndOfEventAction(const G4Event* evt)
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {
   G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();

   G4int n_trajectories = 0;

   if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
   for(G4int i=0; i<n_trajectories; i++)
      { G4Trajectory* trj = 
                      (G4Trajectory *)((*(evt->GetTrajectoryContainer()))[i]);
        if (fDrawFlag == "all") trj->DrawTrajectory(50);
        else if ((fDrawFlag == "charged")&&(trj->GetCharge() != 0.))
                               trj->DrawTrajectory(50);
      }
  }

  if (fVerboselevel>0)
     G4cout << "<<< Event  " << evt->GetEventID() << " ended." << G4endl;

  if (fRunaction->GetRndmFreq() == 2)
    {
     CLHEP::HepRandom::saveEngineStatus("endOfEvent.rndm");
     G4int evtNb = evt->GetEventID();
     if (evtNb%fPrintModulo == 0)
       {
        G4cout << "\n---> End of Event: " << evtNb << G4endl;
        CLHEP::HepRandom::showEngineStatus();
       }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int F04EventAction::GetEventNo()
{
  G4int evno = fpEventManager->GetConstCurrentEvent()->GetEventID();
  return evno ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04EventAction::SetEventVerbose(G4int level)
{
  fVerboselevel = level ;
}
