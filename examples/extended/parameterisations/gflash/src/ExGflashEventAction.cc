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
// Created by Joanna Weng 26.11.2004


#include "ExGflashEventAction.hh"
#include "ExGflashHit.hh"
#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Event.hh"
//std
#include <iostream>
#include <algorithm>
//Gflash
using namespace std;

#include "G4Timer.hh"
extern G4Timer Timer;
extern G4Timer Timerintern;



ExGflashEventAction::ExGflashEventAction():
nevent(0),dtime(0.0),calorimeterCollectionId(-1)
{
}

ExGflashEventAction::~ExGflashEventAction()
{
	G4cout << "Internal Real Elapsed Time /event is: "<< dtime /nevent<< G4endl;
}


void ExGflashEventAction::BeginOfEventAction(const G4Event *evt){
	Timerintern.Start();
	G4cout<<" ------ Start ExGflashEventAction ----- "<<G4endl;
	nevent=evt->GetEventID();
	G4cout<<" Start generating event Nr "<<nevent<<G4endl<<G4endl; 	
}

void ExGflashEventAction::EndOfEventAction(const G4Event *evt)
{  
	Timerintern.Stop();
	G4cout << G4endl;
	G4cout << "******************************************";
	G4cout << G4endl;
	G4cout << "Internal Real Elapsed Time is: "<< Timerintern.GetRealElapsed();
	G4cout << G4endl;
	G4cout << "Internal System Elapsed Time: " << Timerintern.GetSystemElapsed();
	G4cout << G4endl;
	G4cout << "Internal GetUserElapsed Time: " << Timerintern.GetUserElapsed();
	G4cout << G4endl;
	G4cout << "******************************************"<< G4endl;
	dtime+=Timerintern.GetRealElapsed();
	G4cout<<" ------ ExGflashEventAction::End of event nr. "<<nevent<<"  -----"<< G4endl;     

	G4SDManager * SDman = G4SDManager::GetSDMpointer();
	G4String colNam;
	calorimeterCollectionId=SDman->GetCollectionID(colNam="ExGflashCollection");
	if (calorimeterCollectionId<0) return;
	G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
	ExGflashHitsCollection* THC = 0;
	G4double totE = 0;
	// Read out of the crysta ECAL
	THC=(ExGflashHitsCollection *)(HCE->GetHC(calorimeterCollectionId));
	if (THC) {
		/// Hits in sensitive Detector
		int n_hit = THC->entries();
		G4cout<<"  " << n_hit<< " hits are stored in ExGflashHitsCollection "<<G4endl;
		G4PrimaryVertex* pvertex=evt->GetPrimaryVertex();   
		///Computing (x,y,z) of vertex of initial particles  
		G4ThreeVector vtx=pvertex->GetPosition();
		G4PrimaryParticle* pparticle=pvertex->GetPrimary();
		// direction of the Shower
		G4ThreeVector mom=pparticle->GetMomentum()/pparticle->GetMomentum().mag();
		
		//@@@ ExGflashEventAction: Magicnumber
		G4double energyincrystal[100];
		for (int i=0;i<100;i++) energyincrystal[i]=0.;
		
		//@@@ ExGflashEventAction: Magicnumber
		/// For all Hits in sensitive detector
		for (int i=0;i<n_hit;i++)
		{
			G4double estep = (*THC)[i]->GetEdep()/GeV;
			if (estep >0.0)
			{
				totE += (*THC)[i]->GetEdep()/GeV;
				G4int num=(*THC)[i]->GetCrystalNum();
				
				energyincrystal[num]+=(*THC)[i]->GetEdep()/GeV;
				//G4cout << num << G4endl;
			//	G4cout << " Crystal Nummer " <<  (*THC)[i]->GetCrystalNum()  << G4endl;
			//	G4cout <<  (*THC)[i]->GetCrystalNum() /10 << "  "<<(*THC)[i]->GetCrystalNum()%10 << G4endl;
				
				G4ThreeVector hitpos=(*THC)[i]->GetPos();	  				
				G4ThreeVector l (hitpos.x(), hitpos.y(), hitpos.z());
				// distance from shower start
				l = l - vtx; 
				// projection on shower axis = longitudinal profile
				G4ThreeVector longitudinal  =  l;	
				// shower profiles (Radial)
				G4ThreeVector radial = vtx.cross(l);
			}
		}
		G4double max=0;
		G4int index=0;
		//Find crystal with maximum energy
		for (int i=0;i<100;i++) 
		{
			//G4cout << i <<"  " << energyincrystal[i] << G4endl;
			if (max <energyincrystal[i])
			{
				max=energyincrystal[i];
				index = i;
			}
		}	
	//G4cout << index <<" NMAX  " << index << G4endl;	

	// 3x3 matrix of crystals around the crystal with the maximum energy deposit
	G4double e3x3 = energyincrystal[index]+energyincrystal[index+1]+energyincrystal[index-1]+
	energyincrystal[index-10]+energyincrystal[index-9]+energyincrystal[index-11]+
	energyincrystal[index+10]+energyincrystal[index+11]+energyincrystal[index+9];

	// 5x5 matrix of crystals around the crystal with the maximum energy deposit	
	G4double e5x5 = energyincrystal[index]+energyincrystal[index+1]+energyincrystal[index-1]+energyincrystal[index+2]+energyincrystal[index-2]+
	energyincrystal[index-10]+energyincrystal[index-9]+energyincrystal[index-11]+energyincrystal[index-8]+energyincrystal[index-12]+
	energyincrystal[index+10]+energyincrystal[index+11]+energyincrystal[index+9]+energyincrystal[index+12]+energyincrystal[index+8];
	
	G4cout << "   e1  " << energyincrystal[index]  << "   e3x3  " << e3x3<<  "   GeV  e5x5  "   <<e5x5  <<G4endl;	
	}
	
	G4cout << " Total energy deposited in the calorimeter: " << totE << " (GeV)" << G4endl;
	G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
	G4int n_trajectories = 0;
	if(trajectoryContainer){ n_trajectories = trajectoryContainer->entries(); }
	G4cout << "    " << n_trajectories  << " trajectories stored in this event." << G4endl;
}












