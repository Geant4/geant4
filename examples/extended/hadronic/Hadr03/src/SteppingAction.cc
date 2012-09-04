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
/// \file hadronic/Hadr03/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// $Id: SteppingAction.cc,v 1.7 2010-10-13 13:42:33 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"
#include "G4HadronicProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(PrimaryGeneratorAction* prim,
                               RunAction* RuAct, HistoManager* Hist)
:fPrimary(prim),fRunAction(RuAct), fHistoManager(Hist)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // count processes
  // 
  const G4StepPoint* endPoint = aStep->GetPostStepPoint();
  const G4VProcess* process   = endPoint->GetProcessDefinedStep();
  fRunAction->CountProcesses(process);
  
  // check that an real interaction occured (eg. not a transportation)     
  G4bool transmit = (endPoint->GetStepStatus() <= fGeomBoundary);  
  if (transmit) return;
                      
  //real processes : sum track length
  //
  G4double stepLength = aStep->GetStepLength();
  fRunAction->SumTrack(stepLength);
  
  //energy-momentum balance initialisation
  //
  const G4StepPoint* prePoint = aStep->GetPreStepPoint();
  G4double Q             = - prePoint->GetKineticEnergy();
  G4ThreeVector Pbalance = - prePoint->GetMomentum();
  
  //initialisation of the nuclear channel identification
  //
  G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition();
  G4String partName = particle->GetParticleName();
  G4String nuclearChannel = partName;
  G4HadronicProcess* hproc = (G4HadronicProcess*) process;  
  G4String targetName = hproc->GetTargetIsotope()->GetName();
  nuclearChannel += " + " + targetName + " --> ";
    
  //scattered primary particle (if any)
  //
  G4int ih = 1;
  if (aStep->GetTrack()->GetTrackStatus() == fAlive) {
    G4double energy = endPoint->GetKineticEnergy();      
    fHistoManager->FillHisto(ih,energy);
    //
    G4ThreeVector momentum = endPoint->GetMomentum();
    Q        += energy;
    Pbalance += momentum;
    //
    nuclearChannel += partName + " + ";
  }  
  
  //secondaries
  //
  const G4TrackVector* secondary = fpSteppingManager->GetSecondary();
  for (size_t lp=0; lp<(*secondary).size(); lp++) {
    particle = (*secondary)[lp]->GetDefinition(); 
    G4String name   = particle->GetParticleName();
    G4String type   = particle->GetParticleType();      
    G4double charge = particle->GetPDGCharge();
    G4double energy = (*secondary)[lp]->GetKineticEnergy();
    fRunAction->ParticleCount(name,energy);
    //energy spectrum
    if (charge > 3.)  ih = 2; 
    else if (particle == G4Gamma::Gamma())       ih = 3;
    else if (particle == G4Neutron::Neutron())   ih = 4;
    else if (particle == G4Proton::Proton())     ih = 5;
    else if (particle == G4Deuteron::Deuteron()) ih = 6;
    else if (particle == G4Alpha::Alpha())       ih = 7;       
    else if (type == "nucleus")                  ih = 8;
    else if (type == "meson")                    ih = 9;
    else if (type == "baryon")                   ih = 10;        
    fHistoManager->FillHisto(ih,energy);
    //energy-momentum balance
    G4ThreeVector momentum = (*secondary)[lp]->GetMomentum();
    Q        += energy;
    Pbalance += momentum;
    //particle flag
    fParticleFlag[particle] = true;
  }
  
  //energy-momentum balance
  G4double Pbal = Pbalance.mag();
  fRunAction->Balance(Pbal);
  ih = 11;
  fHistoManager->FillHisto(ih,Q);
  ih = 12;
  fHistoManager->FillHisto(ih,Pbal);  
  
  // nuclear channel
 std::map<G4ParticleDefinition*,G4bool>::iterator ip;               
 for (ip = fParticleFlag.begin(); ip != fParticleFlag.end(); ip++) { 
    G4String name = ip->first->GetParticleName();
    if (ip != fParticleFlag.begin()) nuclearChannel += " + ";
    nuclearChannel += name;
 }
 
  ///G4cout << "\n nuclear channel: " << nuclearChannel << G4endl;
  fRunAction->CountNuclearChannel(nuclearChannel, Q);
    
  fParticleFlag.clear();
              
  // kill event after first interaction
  //
  G4RunManager::GetRunManager()->AbortEvent();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


