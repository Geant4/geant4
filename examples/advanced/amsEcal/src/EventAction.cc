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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "Run.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* det,PrimaryGeneratorAction* prim)
:detector(det), primary(prim)
{
  nbOfModules = detector->GetNbModules();	 	
  nbOfLayers  = detector->GetNbLayers();
  kLayerMax = nbOfModules*nbOfLayers + 1;
  
  EtotCalor = EvisCalor = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
  EtotLayer.resize(kLayerMax);
  EvisLayer.resize(kLayerMax);			
  for (G4int k=0; k<kLayerMax; k++) {
    EtotLayer[k] = EvisLayer[k] = 0.0;
  }
  EtotCalor = EvisCalor = 0.;
  EvisFiber.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::SumDeStep(G4int iModule, G4int iLayer, G4int iFiber,
                            G4double deStep )
{
  if (iModule > 0) EtotCalor += deStep;
  		
  G4int kLayer = 0; G4int kFiber = 0;
  if (iLayer > 0) {
	kLayer = (iModule-1)*nbOfLayers + iLayer;
	EtotLayer[kLayer] += deStep;
  }
  
  if (iFiber > 0) {
	EvisLayer[kLayer] += deStep;
	EvisCalor += deStep;
    kFiber = 1000*kLayer + iFiber;
	EvisFiber[kFiber] += deStep;	  	
  }	  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{
  //pass informations to Run
  //
  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  	
  for (G4int k=0; k<kLayerMax; k++) {
     run->SumEvents_1(k,EtotLayer[k],EvisLayer[k]);   
  }
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();    
  analysisManager->FillH1(1,EtotCalor);
  analysisManager->FillH1(2,EvisCalor);
  
  G4double Ebeam = primary->GetParticleGun()->GetParticleEnergy();
  G4double Eleak = Ebeam - EtotCalor;
  run->SumEvents_2(EtotCalor,EvisCalor,Eleak);
  
  
  std::map<G4int,G4double>::iterator it;         
  for (it = EvisFiber.begin(); it != EvisFiber.end(); it++) {
     G4int kFiber = it->first;
	 G4int iFiber = kFiber%1000;
     G4double Evis = it->second;
	 analysisManager->FillH1(5,iFiber+0.5,Evis);
  }
    
  //write fired fibers on a file
  //
  //// WriteFibers(evt); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
        
#include <fstream>

void EventAction::WriteFibers(const G4Event* evt)
{
  // event is appended on a file
  //
  G4String name = G4AnalysisManager::Instance()->GetFileName();
  G4String fileName = name + ".fibers.ascii";
  
  std::ofstream File(fileName, std::ios::app);
  std::ios::fmtflags mode = File.flags();  
  File.setf( std::ios::scientific, std::ios::floatfield );
  G4int prec = File.precision(3);
    
  //write event number  
  //
  File << evt->GetEventID() << G4endl;
  
  //gun particle informations
  //
  G4ParticleGun* gun = primary->GetParticleGun();
  G4double ekin = gun->GetParticleEnergy();
  G4ThreeVector direction = gun->GetParticleMomentumDirection();
  G4ThreeVector position  = gun->GetParticlePosition();
  File << ekin << " " << direction << " " << position << G4endl;  
  
  //write fibers
  //
  File << EvisFiber.size() << G4endl;
  //
  std::map<G4int,G4double>::iterator it;         
  for (it = EvisFiber.begin(); it != EvisFiber.end(); it++) {
     G4int kFiber = it->first;
     G4double Evis = it->second;
     File << " " << std::setw(7) << kFiber << " "<< std::setw(10) << Evis
            << G4endl;
  }
           
  File << G4endl;
    
  // restaure default formats
  File.setf(mode,std::ios::floatfield);
  File.precision(prec);         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

