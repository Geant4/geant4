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
// $Id: Tst50EventAction.cc,v 1.8 2003-01-28 08:57:50 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4UnitsTable.hh" 
#include "Tst50EventAction.hh"
#ifdef G4ANALYSIS_USE
#include "Tst50AnalysisManager.hh"
#endif
#include "Tst50PrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "Tst50TrackerHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
Tst50EventAction::Tst50EventAction(Tst50PrimaryGeneratorAction* Primary,G4bool RY,G4String file):
  hit_CollID(-1),p_Primary(Primary),RadiationY(RY),filename(file)
{
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
Tst50EventAction::~Tst50EventAction()
{
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void Tst50EventAction::BeginOfEventAction(const G4Event*)
{
  // G4cout<<"evento:"<<GetEventno()<<G4endl;
  energy=0;

if (hit_CollID==-1)
{
      G4SDManager * SDman = G4SDManager::GetSDMpointer();
      hit_CollID = SDman->GetCollectionID("trackerCollection");
      //the pointer points to the ID number of the sensitive detector
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void Tst50EventAction::EndOfEventAction(const G4Event* evt)
{
 G4double initialEnergy=p_Primary->GetInitialEnergy();
  
  // extracted from hits, compute the total energy deposit 

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  
  Tst50TrackerHitsCollection* hit_HC = 0;
  G4int n_hit = 0;
  G4double energyDep=0.;
  
  if (HCE) hit_HC = (Tst50TrackerHitsCollection*)(HCE->GetHC(hit_CollID));
  if(hit_HC)
    { 
      n_hit = hit_HC->entries();
      for (G4int i=0;i<n_hit;i++)
	{ 
	   energyDep= (((*hit_HC)[i]->GetEdep()))/MeV; 
	
/********         
#ifdef G4ANALYSIS_USE
	    Tst50AnalysisManager* analysis = Tst50AnalysisManager::getInstance();
	    if(energyDep!=0)
	        analysis->energy_deposit(energyDep);
#endif
*/
	}
    }
  /*
G4cout<<"energia iniziale in MeV:"<<initialEnergy/MeV<<G4endl;
G4cout<<"energia in MeV:"<<energy/MeV<<G4endl;
  G4double radiation=(energy/initialEnergy);
  G4cout<<"Radiation yield:"<< radiation<<G4endl;
 
  */
  if (RadiationY)
    {
G4std::ofstream pmtfile(filename, G4std::ios::app);
 G4double radiation=(energy/MeV)/(initialEnergy/MeV);
 G4cout<<"Radiation yield:"<< radiation<<G4endl;
if(pmtfile.is_open()){
  pmtfile<<'\t'<<radiation<<'\t'<<'\t'<<initialEnergy/MeV<<G4endl;
     }
    }
 // get number of stored trajectories
  //
  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  
  // periodic printing
  //
      


  // extract the trajectories and draw them
  //
  if (G4VVisManager::GetConcreteInstance())
    {
     for (G4int i=0; i<n_trajectories; i++) 
        { G4Trajectory* trj = (G4Trajectory*)
	                            ((*(evt->GetTrajectoryContainer()))[i]);
          trj->DrawTrajectory(50);
        }
    }
}
G4int Tst50EventAction::GetEventno()
{
  G4int evno = fpEventManager->GetConstCurrentEvent()->GetEventID() ;
  return evno ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double Tst50EventAction::RadiationYield(G4double energyLost)
{
  energy += energyLost;
}


