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
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelEventAction.cc                           *     
// * -------                                                            *
// *                                                                    *
// * Version:           0.4                                             *
// * Date:              06/11/00                                        *
// * Author:            R Nartallo                                      *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
// 
// CHANGE HISTORY
// --------------
//
// 05.12.2001 R.Nartallo
// - Implement "Update" function as in LowEnTest
//
// 06.11.2000 R.Nartallo
// - First implementation of xray_telescope event action
// - Based on Chandra and XMM models 
//
//
// **********************************************************************

#include "G4ios.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"

#include "XrayTelEventAction.hh"
#include "XrayTelEventActionMessenger.hh"


XrayTelEventAction::XrayTelEventAction()
  : drawFlag("all"),eventMessenger(0),drawEvent(false)
{
  eventMessenger = new XrayTelEventActionMessenger(this);
}


XrayTelEventAction::~XrayTelEventAction()
{
  delete eventMessenger;
}


void XrayTelEventAction::BeginOfEventAction(const G4Event* Ev)
{
  drawEvent = false;
}


void XrayTelEventAction::EndOfEventAction(const G4Event* Ev)
{
  if (drawEvent){
    const G4Event* evt = fpEventManager->GetConstCurrentEvent();

    G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if ( trajectoryContainer ){ 
      n_trajectories = trajectoryContainer->entries(); 
    }
    
    if ( G4VVisManager::GetConcreteInstance() ) {
      for ( G4int i=0; i<n_trajectories; i++ ) {
	G4Trajectory* trj = (G4Trajectory*)(*(evt->GetTrajectoryContainer()))[i];
	if ( drawFlag == "all" ) trj->DrawTrajectory(50);
	else if ( (drawFlag == "charged")&&(trj->GetCharge() > 0.) )
	  trj->DrawTrajectory(50); 
	  trj->ShowTrajectory(); 
      }
    }
  }
}


void XrayTelEventAction::Update()
{
  // Update the flag identifying the event to be visualised
  drawEvent = true;
}












































