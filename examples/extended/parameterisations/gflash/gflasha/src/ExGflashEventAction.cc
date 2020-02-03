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
#include "ExGflashHistoManager.hh"
#include "ExGflashDetectorConstruction.hh"
#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
//std
#include <iostream>
#include <algorithm>
#include <vector>
typedef std::vector<G4double> MyVector;

//Gflash
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExGflashEventAction::ExGflashEventAction(ExGflashDetectorConstruction* det)
 : G4UserEventAction(),fCalorimeterCollectionId(-1),fDet(det)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExGflashEventAction::~ExGflashEventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExGflashEventAction::BeginOfEventAction(const G4Event * /* evt */)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExGflashEventAction::EndOfEventAction(const G4Event *evt)
{  

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
      //      G4cout<<"  " << n_hit<< " hits are stored in ExGflashHitsCollection "<<G4endl;
      G4PrimaryVertex* pvertex=evt->GetPrimaryVertex();   
      //Computing (x,y,z) of vertex of initial particles  
      G4ThreeVector vtx=pvertex->GetPosition();
      G4PrimaryParticle* pparticle=pvertex->GetPrimary();
      // direction of the Shower
      G4ThreeVector mom=pparticle->GetMomentumDirection();

      G4double Ekin=pparticle->GetKineticEnergy();
      G4double mass=pparticle->GetParticleDefinition()->GetPDGMass();
      G4double Etot = Ekin/MeV + mass/MeV;

      G4int nLbin = fDet->GetnLtot();
      G4int nRbin = fDet->GetnRtot();
      G4double dLradl = fDet->GetdLradl();
      G4double dRradl = fDet->GetdRradl();

      G4double SDRadl = fDet->GetSDRadLen(); // SD matrial
      // init to to 0.0
      MyVector dEdL(nLbin,0.0);
      MyVector dEdR(nRbin,0.0);

      G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();

      fAnalysisManager->FillH1(1,n_hit + 0.5);
      /// For all Hits in sensitive detector
      for (int i=0;i<n_hit;i++)
        {
          G4double estep = (*THC)[i]->GetEdep();
          fAnalysisManager->FillH1(2,estep/MeV);
          estep /= MeV;
        
          if (estep >0.0)
            {
              totE += estep; 
              
              G4ThreeVector hitpos=(*THC)[i]->GetPos();
              // in shower coordinate system
              // from shower start
              G4ThreeVector l = hitpos - vtx;
              //longitudinal profile
              G4ThreeVector longitudinal  =  l.dot(mom) * mom;
              // shower profiles (Radial)
              G4ThreeVector radial = l - longitudinal;

              G4int SlideNb = G4int((longitudinal.mag()/SDRadl) / dLradl);
              G4int RingNb  = G4int((radial.mag()/SDRadl) / dRradl);

              if ( SlideNb >=0 && SlideNb < nLbin) dEdL[SlideNb] += estep;
              if ( RingNb >=0 && RingNb < nLbin)dEdR[RingNb] += estep;
            }
        }

      G4double Lnorm = 100. / dLradl / Etot;
      G4double Lsum = 0.0;
      for (int i=0;i<nLbin;i++)
        { // Slide
          //      fAnalysisManager->FillH1(3,(i +0.5) * dLradl,dEdL[i] * Lnorm);
          fAnalysisManager->FillP1(0,(i +0.5) * dLradl,dEdL[i] * Lnorm);
          Lsum += dEdL[i];
          fAnalysisManager->FillP1(2,(i +0.5) * dLradl,Lsum * Lnorm);
        }
      G4double Rnorm = 100. / dRradl / Etot;
      G4double Rsum = 0.0;
      for (int i=0;i<nRbin;i++)
        { // Ring
          //      fAnalysisManager->FillH1(4,(i +0.5) * dRradl,dEdR[i] * Rnorm);
          fAnalysisManager->FillP1(1,(i +0.5) * dRradl,dEdR[i] * Rnorm);
          Rsum += dEdR[i];
          fAnalysisManager->FillP1(3,(i +0.5) * dRradl,Rsum * Rnorm);
        }

      fAnalysisManager->FillH1(0,totE/Etot * 100.);

    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
