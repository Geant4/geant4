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

#include "exGPSEventAction.hh"

#ifdef G4ANALYSIS_USE
#include <AIDA/AIDA.h>
#include "exGPSAnalysisManager.hh"
#endif

#include "exGPSEventActionMessenger.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
exGPSEventAction::exGPSEventAction()
:drawFlag("all"), printModulo(1000), eventMessenger(NULL) 
{
  eventMessenger = new exGPSEventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exGPSEventAction::~exGPSEventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exGPSEventAction::BeginOfEventAction(const G4Event* evt)
{
  
  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0)
    { 
      G4cout << "\n---> Begin of event: " << evtNb << G4endl;
      //  HepRandom::showEngineStatus();
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exGPSEventAction::EndOfEventAction(const G4Event* evt)
{
  //  G4int evtNb = evt->GetEventID();

#ifdef G4ANALYSIS_USE 
  G4int nVertex = evt->GetNumberOfPrimaryVertex();
  for ( G4int j = 0; j < nVertex; j++) { 
    G4int nPart =  evt->GetPrimaryVertex(j)->GetNumberOfParticle(); 
    for ( G4int i = 0; i < nPart; i++) {
      G4ThreeVector position=evt->GetPrimaryVertex(j)->GetPosition();
      G4ThreeVector direction=evt->GetPrimaryVertex(j)->GetPrimary(i)->GetMomentum();
      G4double P=direction.mag();
      G4double E=evt->GetPrimaryVertex(j)->GetPrimary(i)->GetMass();
      G4double E0=evt->GetPrimaryVertex(j)->GetPrimary(i)->GetG4code()->GetPDGMass();
      E=std::sqrt(P*P+E0*E0);   
      E -= E0;
      G4String pname = evt->GetPrimaryVertex(j)->GetPrimary(i)->GetG4code()->GetParticleName();
      //
      direction=direction*(1./direction.mag());                
      direction = -direction;  // reverse the direction
      G4double x,y,z,w;
      G4double Phi=direction.phi();
      if (Phi <0) Phi=Phi+twopi;
      G4double Theta=direction.theta();
      x=position.x();
      y=position.y();
      z=position.z();
      w = evt->GetPrimaryVertex(j)->GetPrimary(i)->GetWeight();
      exGPSAnalysisManager::getInstance()->Fill(pname,E,x,y,z,Theta,Phi,w);
    }
  }   
#endif
  
  // extract the trajectories and draw them
  if (G4VVisManager::GetConcreteInstance())
    {
     G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
     G4int n_trajectories = 0;
     if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
     
     for (G4int i=0; i<n_trajectories; i++) 
       { G4Trajectory* trj = (G4Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
       if (drawFlag == "all") trj->DrawTrajectory(0);
       else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
	 trj->DrawTrajectory(0);
       else if ((drawFlag == "neutral")&&(trj->GetCharge() == 0.))
	 trj->DrawTrajectory(0);				   
       }
   }             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....






