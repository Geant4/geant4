// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN03EventAction.cc,v 1.10 2000-06-05 09:54:37 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "ExN03EventAction.hh"

#include "ExN03CalorHit.hh"
#include "ExN03EventActionMessenger.hh"

#include "g4rw/tvordvec.h"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ExN03EventAction::ExN03EventAction()
:calorimeterCollID(-1),drawFlag("all"),eventMessenger(NULL),
 printModulo(1)
{
  eventMessenger = new ExN03EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ExN03EventAction::~ExN03EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExN03EventAction::BeginOfEventAction(const G4Event* evt)
{
  
 G4int evtNb = evt->GetEventID();
 if (evtNb%printModulo == 0)
   { 
    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
    HepRandom::showEngineStatus();
   }
    
 if (calorimeterCollID==-1)
   {
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    calorimeterCollID = SDman->GetCollectionID("CalCollection");
   } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExN03EventAction::EndOfEventAction(const G4Event* evt)
{
  G4int evtNb = evt->GetEventID();

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  ExN03CalorHitsCollection* CHC = NULL;
  G4int n_hit = 0;
  G4double totEAbs=0, totLAbs=0, totEGap=0, totLGap=0;
    
  if (HCE) CHC = (ExN03CalorHitsCollection*)(HCE->GetHC(calorimeterCollID));

  if (CHC)
    {
     n_hit = CHC->entries();
     for (G4int i=0;i<n_hit;i++)
	{
	  totEAbs += (*CHC)[i]->GetEdepAbs(); 
	  totLAbs += (*CHC)[i]->GetTrakAbs();
	  totEGap += (*CHC)[i]->GetEdepGap(); 
	  totLGap += (*CHC)[i]->GetTrakGap();
	}
     }	
	  
  if (evtNb%printModulo == 0) {
    G4cout << "---> End of event: " << evtNb << G4endl;	

    G4cout
       << "   Absorber: total energy: " << G4std::setw(7) << G4BestUnit(totEAbs,"Energy")
       << "       total track length: " << G4std::setw(7) << G4BestUnit(totLAbs,"Length")
       << G4endl
       << "        Gap: total energy: " << G4std::setw(7) << G4BestUnit(totEGap,"Energy")
       << "       total track length: " << G4std::setw(7) << G4BestUnit(totLGap,"Length")
       << G4endl;
	  
    G4cout << "\n     " << n_hit
	   << " hits are stored in ExN03CalorHitsCollection." << G4endl;  	     
  }          
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
