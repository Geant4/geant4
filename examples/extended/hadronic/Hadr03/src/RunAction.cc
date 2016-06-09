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
// $Id$
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
  : fDetector(det), fPrimary(prim)
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

  fTotalCount = fGammaCount = 0;
  fSumTrack = fSumTrack2 = 0.;
  for (G4int i=0; i<3; i++) { fPbalance[i] = 0. ; } 
  for (G4int i=0; i<3; i++) { fNbGamma[i] = 0 ; }
       
  //histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    analysisManager->OpenFile();
  }     
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

void RunAction::CountNuclearChannel(G4String name, G4double Q)
{
  fNuclChannelCount[name]++;
  fNuclChannelQ[name] += Q;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::Balance(G4double Pbal)
{ 
  fPbalance[0] += Pbal;
  //update min max   
  if (fTotalCount == 1) fPbalance[1] = fPbalance[2] = Pbal;  
  if (Pbal < fPbalance[1]) fPbalance[1] = Pbal;
  if (Pbal > fPbalance[2]) fPbalance[2] = Pbal;    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::CountGamma(G4int nGamma)
{ 
  fGammaCount++;
  fNbGamma[0] += nGamma;
  //update min max   
  if (fGammaCount == 1) fNbGamma[1] = fNbGamma[2] = nGamma;  
  if (nGamma < fNbGamma[1]) fNbGamma[1] = nGamma;
  if (nGamma > fNbGamma[2]) fNbGamma[2] = nGamma;    
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
  G4cout << "\n Process calls frequency --->";
  std::map<const G4VProcess*,G4int>::iterator it;    
  for (it = fProcCounter.begin(); it != fProcCounter.end(); it++) {
     G4String procName = it->first->GetProcessName();
     G4int    count    = it->second;
     G4cout << "\t" << procName << "= " << count;
     if (procName == "Transportation") survive = count;
  }
      
  if (survive > 0) {
    G4cout << "\n\n Nb of incident particles surviving after "
           << G4BestUnit(fDetector->GetSize(),"Length") << " of "
           << material->GetName() << " : " << survive << G4endl;
  }
  
  if (fTotalCount == 0) fTotalCount = 1;   //force printing anyway
  
  //compute mean free path and related quantities
  //
  G4double MeanFreePath = fSumTrack /fTotalCount;     
  G4double MeanTrack2   = fSumTrack2/fTotalCount;     
  G4double rms = std::sqrt(std::fabs(MeanTrack2 - MeanFreePath*MeanFreePath));
  G4double CrossSection = 0.0;
  if(MeanFreePath > 0.0) { CrossSection = 1./MeanFreePath; }
  G4double massicMFP = MeanFreePath*density;
  G4double massicCS  = 0.0;
  if(massicMFP > 0.0) { massicCS = 1./massicMFP; }
   
  G4cout << "\n\n MeanFreePath:\t"   << G4BestUnit(MeanFreePath,"Length")
         << " +- "                   << G4BestUnit( rms,"Length")
         << "\tmassic: "             << G4BestUnit(massicMFP, "Mass/Surface")
         << "\n CrossSection:\t"     << CrossSection*cm << " cm^-1 "
         << "\t\tmassic: "         << G4BestUnit(massicCS, "Surface/Mass")
         << G4endl;
         
  //check cross section from G4HadronicProcessStore
  //
  G4cout << "\n Verification : "
         << "crossSections from G4HadronicProcessStore:";
  
  G4HadronicProcessStore* store = G4HadronicProcessStore::Instance();  
  G4double sumc = 0.0;  
  for (it = fProcCounter.begin(); it != fProcCounter.end(); it++) {
    const G4VProcess* process = it->first;
    G4double xs =
    store->GetCrossSectionPerVolume(particle,energy,process,material);                   
    G4double massSigma = xs/density;                                                  
    sumc += massSigma;
    G4String procName = process->GetProcessName();    
    G4cout << "\n    " << procName << "= " 
           << G4BestUnit(massSigma, "Surface/Mass");
  }             
  G4cout << "\n    total = " 
         << G4BestUnit(sumc, "Surface/Mass") << G4endl;
         
 //nuclear channel count
 //
 G4cout << "\n   List of nuclear reactions: \n" << G4endl;
     
 std::map<G4String,G4int>::iterator ic;               
 for (ic = fNuclChannelCount.begin(); ic != fNuclChannelCount.end(); ic++) { 
    G4String name = ic->first;
    G4int count   = ic->second;
    G4double Q = fNuclChannelQ[name]/count;
         
    G4cout << "  " << std::setw(50) << name << ": " << std::setw(7) << count
           << "   Q = " << std::setw(wid) << G4BestUnit(Q, "Energy")
           << G4endl;           
 } 
 
 //Gamma count
 //
 if (fGammaCount > 0) {       
   G4cout << "\n" << std::setw(58) << "Number of gamma: N = " 
           << fNbGamma[1] << " --> " << fNbGamma[2] << G4endl;
 }
       
 //particles count
 //
 G4cout << "\n   List of generated particles: \n" << G4endl;
     
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
 
 //energy momentum balance
 //
 if (fTotalCount > 1) {
    G4double Pbmean = fPbalance[0]/fTotalCount;           
    G4cout << "\n   Momentum balance: Pmean = " 
           << std::setw(wid) << G4BestUnit(Pbmean, "Energy")
           << "\t( "  << G4BestUnit(fPbalance[1], "Energy")
           << " --> " << G4BestUnit(fPbalance[2], "Energy")
           << ") \n" << G4endl;
 }
                  
  //restore default format         
  G4cout.precision(dfprec);
           
  // remove all contents in fProcCounter 
  fProcCounter.clear();
  // remove all contents in fNuclChannel 
  fNuclChannelCount.clear();
  fNuclChannelQ.clear();    
  // remove all contents in fParticleCount
  fParticleCount.clear(); 
  fEmean.clear();  fEmin.clear(); fEmax.clear();  
  
  //save histograms      
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
  if ( analysisManager->IsActive() ) {
    analysisManager->Write();
    analysisManager->CloseFile();
  }
      
  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
