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
// $Id: TrackingAction.cc,v 1.2 2010-10-11 14:31:39 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "TrackingAction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingMessenger.hh"

#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(HistoManager* histo,
                               RunAction* RA, EventAction* EA)
:histoManager(histo),run(RA),event(EA)
{
  fullChain = false;
  trackMessenger = new TrackingMessenger(this);   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::~TrackingAction()
{
  delete trackMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  G4ParticleDefinition* particle = track->GetDefinition();
  G4String name   = particle->GetParticleName();
  charge = particle->GetPDGCharge();
  mass   = particle->GetPDGMass();  
    
  G4double Ekin = track->GetKineticEnergy();
  G4int ID      = track->GetTrackID();
  
  G4bool condition = false;  

  //count particles
  //
  run->ParticleCount(name, Ekin);
  
  //energy spectrum
  //
  G4int ih = 0;
  if (particle == G4Electron::Electron()||
      particle == G4Positron::Positron())  ih = 1;
  else if (particle == G4NeutrinoE::NeutrinoE()||
           particle == G4AntiNeutrinoE::AntiNeutrinoE()) ih = 2;
  else if (particle == G4Gamma::Gamma()) ih = 3;
  else if (particle == G4Alpha::Alpha()) ih = 4;
  else if (charge > 2.) ih = 5;
  if (ih) histoManager->FillHisto(ih, Ekin);
  
  //fullChain: stop ion and print decay chain
  //
  if (charge > 2.) {
    G4Track* tr = (G4Track*) track;
    if (fullChain) tr->SetTrackStatus(fStopButAlive);
    if (ID == 1) event->AddDecayChain(name);
      else       event->AddDecayChain(" ---> " + name); 
  }
  
  //example of saving random number seed of this event, under condition
  //
  ////condition = (ih == 3);
  if (condition) G4RunManager::GetRunManager()->rndmSaveThisEvent();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
  //keep only ions
  //
  if (charge < 3. ) return;

  //get time
  //   
  G4double time = track->GetGlobalTime();
    
  //energy and momentum balance (from secondaries)
  //
  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  size_t nbtrk = (*secondaries).size();
  if (nbtrk) {
    //there are secondaries --> it is a decay
    //
    //force 'single' decay
    G4int ID = track->GetTrackID();
    if ((!fullChain)&&(ID > 1)) G4RunManager::GetRunManager()->AbortEvent();
    //
    //balance    
    G4double Ebalance = - track->GetTotalEnergy();
    G4ThreeVector Pbalance = - track->GetMomentum();
    for (size_t itr=0; itr<nbtrk; itr++) {
       G4Track* trk = (*secondaries)[itr];
       Ebalance += trk->GetTotalEnergy();
       //exclude gamma desexcitation from momentum balance
       if (trk->GetDefinition() != G4Gamma::Gamma())	 
         Pbalance += trk->GetMomentum();	         
    }
    G4double Pbal = Pbalance.mag();  
    run->Balance(Ebalance,Pbal);  
    histoManager->FillHisto(6,Ebalance);
    histoManager->FillHisto(7,Pbal);
  }
  
  //no secondaries --> end of chain    
  //  
  if (!nbtrk) {
    run->EventTiming(time);		//time of life
    histoManager->FillHisto(8,time);
    histoManager->FillHisto(9,time);	//activity                    
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

