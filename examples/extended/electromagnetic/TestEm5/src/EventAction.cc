//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: EventAction.cc,v 1.2 2004/02/19 18:18:52 maire Exp $
// GEANT4 tag $Name: geant4-06-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "RunAction.hh"
#include "HistoManager.hh"
#include "EventMessenger.hh"

#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"

#ifdef G4ANALYSIS_USE
 #include "AIDA/IHistogram1D.h"
#endif
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* RA, HistoManager* histo)
:runaction(RA),histoManager(histo),
 drawFlag("none"),printModulo(10000)
{
  eventMessenger = new EventMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{
 G4int evtNb = evt->GetEventID();

 //printing survey
 if (evtNb%printModulo == 0) 
    G4cout << "\n---> Begin of Event: " << evtNb << G4endl;

 // initialisation per event
 EnergyDeposit  = 0.;
 TrakLenCharged = TrakLenNeutral = 0.; 
 nbStepsCharged = nbStepsNeutral = 0;
 TransmitFlag   = ReflectFlag    = 0;    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
 runaction->AddEnergy(EnergyDeposit);
 runaction->AddTrakLenCharg(TrakLenCharged);
 runaction->AddTrakLenNeutr(TrakLenNeutral);
  
 runaction->CountStepsCharg(nbStepsCharged);
 runaction->CountStepsNeutr(nbStepsNeutral);
 
 runaction->CountTransmit (TransmitFlag);
 runaction->CountReflect  (ReflectFlag);
      
#ifdef G4ANALYSIS_USE
 if (histoManager->GetHisto(1)) {
    G4double unit = histoManager->GetHistoUnit(1);
    histoManager->GetHisto(1)->fill(EnergyDeposit/unit);
 }   
#endif

 if (G4VVisManager::GetConcreteInstance())
  {
   G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
   G4int n_trajectories = 0;
   if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();  
   for (G4int i=0; i<n_trajectories; i++) 
      { G4Trajectory* trj = (G4Trajectory*)
                                      ((*(evt->GetTrajectoryContainer()))[i]);
        if (drawFlag == "all") trj->DrawTrajectory(1000);
        else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
                               trj->DrawTrajectory(1000); 
      }
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

