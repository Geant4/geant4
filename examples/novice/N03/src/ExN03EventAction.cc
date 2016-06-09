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
// $Id: ExN03EventAction.cc,v 1.27 2006/06/29 17:49:00 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN03EventAction.hh"

#include "ExN03RunAction.hh"
#include "ExN03EventActionMessenger.hh"

#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03EventAction::ExN03EventAction(ExN03RunAction* run)
:runAct(run),drawFlag("all"),printModulo(1),eventMessenger(0)
{
  eventMessenger = new ExN03EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03EventAction::~ExN03EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03EventAction::BeginOfEventAction(const G4Event* evt)
{  
 G4int evtNb = evt->GetEventID();
 if (evtNb%printModulo == 0) { 
   G4cout << "\n---> Begin of event: " << evtNb << G4endl;
   CLHEP::HepRandom::showEngineStatus();
 }
 
 // initialisation per event
 EnergyAbs = EnergyGap = 0.;
 TrackLAbs = TrackLGap = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03EventAction::EndOfEventAction(const G4Event* evt)
{
  //accumulates statistic
  //
  runAct->fillPerEvent(EnergyAbs, EnergyGap, TrackLAbs, TrackLGap);
  
  //print per event (modulo n)
  //
  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0) {
    G4cout << "---> End of event: " << evtNb << G4endl;	

    G4cout
       << "   Absorber: total energy: " << std::setw(7)
                                        << G4BestUnit(EnergyAbs,"Energy")
       << "       total track length: " << std::setw(7)
                                        << G4BestUnit(TrackLAbs,"Length")
       << G4endl
       << "        Gap: total energy: " << std::setw(7)
                                        << G4BestUnit(EnergyGap,"Energy")
       << "       total track length: " << std::setw(7)
                                        << G4BestUnit(TrackLGap,"Length")
       << G4endl;
	  
  }
  
  // extract the trajectories and draw them

  // You can get a default drawing without this code by using, e.g.,
  // /vis/scene/add/trajectories 1000
  // The code here adds sophistication under control of drawFlag.

  // See comments in G4VTrajectory::DrawTrajectory for the
  // interpretation of the argument, 1000.
  
  G4VVisManager* pVisManager = G4VVisManager::GetConcreteInstance();
  if (pVisManager)
    {
     G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
     G4int n_trajectories = 0;
     if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

     for (G4int i=0; i<n_trajectories; i++) 
        { G4VTrajectory* trj = ((*(evt->GetTrajectoryContainer()))[i]);
          if (drawFlag == "all") pVisManager->Draw(*trj,100);
          else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
                                  pVisManager->Draw(*trj,100);
          else if ((drawFlag == "neutral")&&(trj->GetCharge() == 0.))
                                  pVisManager->Draw(*trj,100);
        }
  }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
