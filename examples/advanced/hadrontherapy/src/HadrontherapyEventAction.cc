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
// $Id: HadrontherapyEventAction.cc,v 1.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------
#include "HadrontherapyEventAction.hh"
#include "HadrontherapyRunAction.hh"
#include "HadrontherapyHit.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "HadrontherapyCalorimeterSD.hh"

// -------------------------------------------------------------------
HadrontherapyEventAction::HadrontherapyEventAction(HadrontherapyRunAction* runAction)
  :calorimeterCollID(-1), eventMessenger(NULL), verboselevel(0)
{
  p_Run = runAction;
}

// --------------------------------------------------------------------
HadrontherapyEventAction::~HadrontherapyEventAction()
{
}

// --------------------------------------------------------------------
void HadrontherapyEventAction::BeginOfEventAction(const G4Event* evt)
{
  event_id = evt->GetEventID();
  if (calorimeterCollID==-1)
    {
      G4SDManager * SDman = G4SDManager::GetSDMpointer();
      calorimeterCollID = SDman->GetCollectionID("CalCollection");
    } 
 
  nstep = 0. ;
  for (G4int slice = 0; slice < 20000; slice ++)      //slice dependence
    {
      energyDep[slice]=0;
    }
}

// ---------------------------------------------------------------------
void HadrontherapyEventAction::EndOfEventAction(const G4Event* evt)
{
  if (calorimeterCollID < 0) return;
  
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  HadrontherapyHitsCollection* CHC = NULL; 
  if(HCE)
    CHC = (HadrontherapyHitsCollection*)(HCE->GetHC(calorimeterCollID));
  
  if(CHC)
    {	
      G4int HitCount = CHC->entries();
      
      for (G4int h=0; h < HitCount; h++)
	{
	  G4int  j=((*CHC)[h])->GetSliceID();  
	  energyDep[j]+=((*CHC)[h]->GetEdep());
	  p_Run->EnergyTotSlice(j,energyDep[j]);	  
	}
    }  
}
G4int HadrontherapyEventAction::Trasporto()
{ 
  G4int var = event_id;
  return var;
}






