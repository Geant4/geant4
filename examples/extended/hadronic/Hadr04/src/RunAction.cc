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
/// \file hadronic/Hadr03/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// $Id: RunAction.cc 70756 2013-06-05 12:20:06Z ihrivnac $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4HadronicProcessStore.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* prim)
  : G4UserRunAction(),
    fDetector(det), fPrimary(prim), fHistoManager(0)
{
 fHistoManager = new HistoManager(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
 delete fHistoManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  CLHEP::HepRandom::showEngineStatus();
       
  //histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    analysisManager->OpenFile();
  }  
  //initialization per run
  //
  fNbStep1 = fNbStep2 = 0;
  fTrackLen1 = fTrackLen2 = 0.;
  fTime1 = fTime2 = 0.;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::ParticleCount(G4String name, G4double Ekin)
{
  fParticleCount[name]++;
  fEmean[name] += Ekin;
  //update min max
  if (fParticleCount[name] == 1) fEmin[name] = fEmax[name] = Ekin;
  if (Ekin < fEmin[name]) fEmin[name] = Ekin;
  if (Ekin > fEmax[name]) fEmax[name] = Ekin;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SumTrackLength(G4int nstep1, G4int nstep2, 
                               G4double trackl1, G4double trackl2,
                               G4double time1, G4double time2)
{
  fNbStep1   += nstep1;  fNbStep2   += nstep2;
  fTrackLen1 += trackl1; fTrackLen2 += trackl2;
  fTime1 += time1; fTime2 += time2;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;

  G4int prec = 5, wid = prec + 2;  
  G4int dfprec = G4cout.precision(prec);
    
  G4Material* material = fDetector->GetMaterial();
  G4double density = material->GetDensity();
  G4int survive = 0;
   
  G4ParticleDefinition* particle = 
                            fPrimary->GetParticleGun()->GetParticleDefinition();
  G4String Particle = particle->GetParticleName();    
  G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
  G4cout << "\n The run consists of " << NbOfEvents << " "<< Particle << " of "
         << G4BestUnit(energy,"Energy") << " through " 
         << G4BestUnit(fDetector->GetSize(),"Length") << " of "
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
  
  //frequency of processes
  //
  G4cout << "\n Process calls frequency --->";
  std::map<const G4VProcess*,G4int>::iterator it;    
  for (it = fProcCounter.begin(); it != fProcCounter.end(); it++) {
     G4String procName = it->first->GetProcessName();
     G4int    count    = it->second;
     G4cout << "\t" << procName << "= " << count;
     if (procName == "Transportation") survive = count;
  }
  G4cout << G4endl;
      
  if (survive > 0) {
    G4cout << "\n Nb of incident particles surviving after "
           << G4BestUnit(0.5*(fDetector->GetSize()),"Length") << " of "
           << material->GetName() << " : " << survive << G4endl;
  }

 // total track length of incident neutron
 //
 G4cout << "\n Parcours of incident neutron:";
  
 G4double meanCollision1  = (G4double)fNbStep1/NbOfEvents;
 G4double meanCollision2  = (G4double)fNbStep2/NbOfEvents; 

 G4cout << "\n   nb of collisions    E>1*eV= " << meanCollision1
        << ";  E<1*eV= " << meanCollision2;
        
 G4double meanTrackLen1  = fTrackLen1/NbOfEvents;
 G4double meanTrackLen2  = fTrackLen2/NbOfEvents; 

 G4cout 
   << "\n   total track length  E>1*eV= " << G4BestUnit(meanTrackLen1,"Length")
   << "  E<1*eV= " << G4BestUnit(meanTrackLen2, "Length");
   
 G4double meanTime1  = fTime1/NbOfEvents;
 G4double meanTime2  = fTime2/NbOfEvents; 

 G4cout 
   << "\n   time of flight      E>1*eV= " << G4BestUnit(meanTime1,"Time")
   << "  E<1*eV= " << G4BestUnit(meanTime2, "Time") << G4endl;
             
 //particles count
 //
 G4cout << "\n List of generated particles:" << G4endl;
     
 std::map<G4String,G4int>::iterator ip;               
 for (ip = fParticleCount.begin(); ip != fParticleCount.end(); ip++) { 
    G4String name = ip->first;
    G4int count   = ip->second;
    G4double eMean = fEmean[name]/count;
    G4double eMin = fEmin[name], eMax = fEmax[name];    
         
    G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
           << "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
           << "\t( "  << G4BestUnit(eMin, "Energy")
           << " --> " << G4BestUnit(eMax, "Energy") 
           << ")" << G4endl;           
 }
                   
  //restore default format         
  G4cout.precision(dfprec);
           
  // remove all contents in fProcCounter 
  fProcCounter.clear();
  // remove all contents in fParticleCount
  fParticleCount.clear(); 
  fEmean.clear();  fEmin.clear(); fEmax.clear();  
  
  //save histograms      
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  ////G4double factor = 1./NbOfEvents;
  ////analysisManager->ScaleH1(3,factor);      
  if ( analysisManager->IsActive() ) {
    analysisManager->Write();
    analysisManager->CloseFile();
  }
      
  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
