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

using namespace std;

#include "Tst34EventAction.hh"
#include "Tst34Hit.hh"
#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4UImanager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Event.hh"
//std
#include <iostream>
#include <algorithm>
//Gflash


#include "G4Timer.hh"
extern G4Timer Timer;
extern G4Timer Timerintern;



Tst34EventAction::Tst34EventAction():
nevent(0),dtime(0.0),calorimeterCollectionId(-1)
{
}

Tst34EventAction::~Tst34EventAction()
{
	cout << "Internal Real Elapsed Time /event is: "<< dtime /nevent<< endl;
}


void Tst34EventAction::BeginOfEventAction(const G4Event *evt){
	Timerintern.Start();
	cout<<" ------ Start Tst34EventAction ----- "<<endl;
	nevent=evt->GetEventID();
	cout<<" Start generating event Nr "<<nevent<<endl<<endl; 	
}

void Tst34EventAction::EndOfEventAction(const G4Event *evt)
{  
	Timerintern.Stop();
	cout << endl;
	cout << "******************************************";
	cout << endl;
	cout << "Internal Real Elapsed Time is: "<< Timerintern.GetRealElapsed();
	cout << endl;
	cout << "Internal System Elapsed Time: " << Timerintern.GetSystemElapsed();
	cout << endl;
	cout << "Internal GetUserElapsed Time: " << Timerintern.GetUserElapsed();
	cout << endl;
	cout << "******************************************"<< endl;
	dtime+=Timerintern.GetRealElapsed();
	cout<<" ------ Tst34EventAction::End of event nr. "<<nevent<<"  -----"<< endl;     

	G4SDManager * SDman = G4SDManager::GetSDMpointer();
	G4String colNam;
	calorimeterCollectionId=SDman->GetCollectionID(colNam="Tst34Collection");
	if (calorimeterCollectionId<0) return;
	G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
	Tst34HitsCollection* THC = 0;
	G4double totE = 0;
	// Read out of the crysta ECAL
	THC=(Tst34HitsCollection *)(HCE->GetHC(calorimeterCollectionId));
	if (THC) {
		/// Hits in sensitive Detector
		int n_hit = THC->entries();
		cout<<"  " << n_hit<< " hits are stored in Tst34HitsCollection "<<endl;
		G4PrimaryVertex* pvertex=evt->GetPrimaryVertex();   
		///Computing (x,y,z) of vertex of initial particles  
		G4ThreeVector vtx=pvertex->GetPosition();
		G4PrimaryParticle* pparticle=pvertex->GetPrimary();
		// direction of the Shower
		G4ThreeVector mom=pparticle->GetMomentum()/pparticle->GetMomentum().mag();
		
		
		/// For all Hits in sensitive detector
		for (int i=0;i<n_hit;i++)
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
	G4cout << " #### Tst34EventAction::Test: Total energy deposited in the calorimeter: " << totE << " (GeV)" << G4endl;

	
	
}












