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
using namespace std;

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


#include "G4Timer.hh"
extern G4Timer Timer;
extern G4Timer Timerintern;



ExGflashEventAction::ExGflashEventAction():
nevent(0),dtime(0.0),calorimeterCollectionId(-1)
{
}

ExGflashEventAction::~ExGflashEventAction()
{
	cout << "Internal Real Elapsed Time /event is: "<< dtime /nevent<< endl;
}


void ExGflashEventAction::BeginOfEventAction(const G4Event *evt){
	Timerintern.Start();
	cout<<" ------ Start ExGflashEventAction ----- "<<endl;
	nevent=evt->GetEventID();
	cout<<" Start generating event Nr "<<nevent<<endl<<endl; 	
}

void ExGflashEventAction::EndOfEventAction(const G4Event *evt)
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
	cout<<" ------ ExGflashEventAction::End of event nr. "<<nevent<<"  -----"<< endl;     

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
		cout<<"  " << n_hit<< " hits are stored in ExGflashHitsCollection "<<endl;
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
				//std::cout << num << std::endl;
				//@@@ExGflashEventAction: was bringen die namespace statt using namespace std ?
			//	std::cout << " Crystal Nummer " <<  (*THC)[i]->GetCrystalNum()  << std::endl;
			//	std::cout <<  (*THC)[i]->GetCrystalNum() /10 << "  "<<(*THC)[i]->GetCrystalNum()%10 << std::endl;
				
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
			//std::cout << i <<"  " << energyincrystal[i] << std::endl;
			if (max <energyincrystal[i])
			{
				max=energyincrystal[i];
				index = i;
			}
		}	
	//std::cout << index <<" NMAX  " << index << std::endl;	

	// 3x3 matrix of crystals around the crystal with the maximum energy deposit
	G4double e3x3 = energyincrystal[index]+energyincrystal[index+1]+energyincrystal[index-1]+
	energyincrystal[index-10]+energyincrystal[index-9]+energyincrystal[index-11]+
	energyincrystal[index+10]+energyincrystal[index+11]+energyincrystal[index+9];

	// 5x5 matrix of crystals around the crystal with the maximum energy deposit	
	G4double e5x5 = energyincrystal[index]+energyincrystal[index+1]+energyincrystal[index-1]+energyincrystal[index+2]+energyincrystal[index-2]+
	energyincrystal[index-10]+energyincrystal[index-9]+energyincrystal[index-11]+energyincrystal[index-8]+energyincrystal[index-12]+
	energyincrystal[index+10]+energyincrystal[index+11]+energyincrystal[index+9]+energyincrystal[index+12]+energyincrystal[index+8];
	
	std::cout << "   e1  " << energyincrystal[index]  << "   e3x3  " << e3x3<<  "   GeV  e5x5  "   <<e5x5  <<std::endl;	
	}
	
	G4cout << " Total energy deposited in the calorimeter: " << totE << " (GeV)" << G4endl;
	G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
	G4int n_trajectories = 0;
	if(trajectoryContainer){ n_trajectories = trajectoryContainer->entries(); }
	G4cout << "    " << n_trajectories  << " trajectories stored in this event." << endl;
}












