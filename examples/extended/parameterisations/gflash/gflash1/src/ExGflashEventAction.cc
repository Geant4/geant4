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
/// \file ExGflashEventAction.cc
/// \brief Implementation of the ExGflashEventAction class
//
// Created by Joanna Weng 26.11.2004


#include "ExGflashEventAction.hh"
#include "ExGflashHit.hh"
#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
//std
#include <iostream>
#include <algorithm>
//Gflash
using namespace std;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExGflashEventAction::ExGflashEventAction()
 : G4UserEventAction(),
   fNevent(0),fDtime(0.0),fCalorimeterCollectionId(-1)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExGflashEventAction::~ExGflashEventAction()
{
  if ( fNevent > 0 ) {
    G4cout << "Internal Real Elapsed Time /event is: "<< fDtime /fNevent<< G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExGflashEventAction::BeginOfEventAction(const G4Event *evt)
{
  fTimerIntern.Start();
  G4cout<<" ------ Start ExGflashEventAction ----- "<<G4endl;
  fNevent=evt->GetEventID();
  G4cout<<" Start generating event Nr "<<fNevent<<G4endl<<G4endl;   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExGflashEventAction::EndOfEventAction(const G4Event *evt)
{  
  fTimerIntern.Stop();
  G4cout << G4endl;
  G4cout << "******************************************";
  G4cout << G4endl;
  G4cout << "Internal Real Elapsed Time is: "<< fTimerIntern.GetRealElapsed();
  G4cout << G4endl;
  G4cout << "Internal System Elapsed Time: " << fTimerIntern.GetSystemElapsed();
  G4cout << G4endl;
  G4cout << "Internal GetUserElapsed Time: " << fTimerIntern.GetUserElapsed();
  G4cout << G4endl;
  G4cout << "******************************************"<< G4endl;
  fDtime+=fTimerIntern.GetRealElapsed();
  G4cout<<" ------ ExGflashEventAction::End of event nr. "<<fNevent<<"  -----"<< G4endl;
  
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4String colNam;
  fCalorimeterCollectionId=SDman->GetCollectionID(colNam="ExGflashCollection");
  if (fCalorimeterCollectionId<0) return;
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  ExGflashHitsCollection* THC = 0;
  G4double totE = 0;
  // Read out of the crysta ECAL
  THC=(ExGflashHitsCollection *)(HCE->GetHC(fCalorimeterCollectionId));
  if (THC)
    {
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
    G4int hitsincrystal[100];
    for (int i=0;i<100;i++) energyincrystal[i]=0.;
    for (int i=0;i<100;i++) hitsincrystal[i]=0.;
    
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
            hitsincrystal[num]++;
            //G4cout << num << G4endl;
            //  G4cout << " Crystal Nummer " <<  (*THC)[i]->GetCrystalNum()  << G4endl;
            //  G4cout <<  (*THC)[i]->GetCrystalNum() /10 <<
            // "  "<<(*THC)[i]->GetCrystalNum()%10 << G4endl;
            
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
    G4double   max = 0;
    G4int    index = 0;
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
    G4double e3x3 =
      energyincrystal[index]+energyincrystal[index+1]+energyincrystal[index-1]+
      energyincrystal[index-10]+energyincrystal[index-9]+energyincrystal[index-11]+
      energyincrystal[index+10]+energyincrystal[index+11]+energyincrystal[index+9];
    
    // 5x5 matrix of crystals around the crystal with the maximum energy deposit  
    G4double e5x5 =
      energyincrystal[index]+energyincrystal[index+1]+energyincrystal[index-1]+
      energyincrystal[index+2]+energyincrystal[index-2]+
      energyincrystal[index-10]+energyincrystal[index-9]+energyincrystal[index-11]+
      energyincrystal[index-8]+energyincrystal[index-12]+
      energyincrystal[index+10]+energyincrystal[index+11]+energyincrystal[index+9]+
      energyincrystal[index+12]+energyincrystal[index+8];

    // 3x3 matrix of crystals around the crystal with the maximum energy deposit
    G4int num3x3 =
      hitsincrystal[index]+hitsincrystal[index+1]+hitsincrystal[index-1]+
      hitsincrystal[index-10]+hitsincrystal[index-9]+hitsincrystal[index-11]+
      hitsincrystal[index+10]+hitsincrystal[index+11]+hitsincrystal[index+9];
    
    // 5x5 matrix of crystals around the crystal with the maximum energy deposit  
    G4int num5x5 =
      hitsincrystal[index]+hitsincrystal[index+1]+hitsincrystal[index-1]+
      hitsincrystal[index+2]+hitsincrystal[index-2]+
      hitsincrystal[index-10]+hitsincrystal[index-9]+hitsincrystal[index-11]+
      hitsincrystal[index-8]+hitsincrystal[index-12]+
      hitsincrystal[index+10]+hitsincrystal[index+11]+hitsincrystal[index+9]+
      hitsincrystal[index+12]+hitsincrystal[index+8];
    
    G4cout << "   e1  " << energyincrystal[index]  
           << "   e3x3  " << e3x3<<  "   GeV  e5x5  "   <<e5x5  <<G4endl;  
    
    G4cout << "   num1  " << hitsincrystal[index]  
           << "   num3x3  " << num3x3<<  "    num5x5  "   <<num5x5  <<G4endl;  
  }
  
  G4cout << " Total energy deposited in the calorimeter: " << totE << " (GeV)" << G4endl;
  G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if(trajectoryContainer){ n_trajectories = trajectoryContainer->entries(); }
  G4cout << "    " << n_trajectories  << " trajectories stored in this event." << G4endl;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......













