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
/// \file radioactivedecay/rdecay01/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
//
// $Id$
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "RunAction.hh"
#include "HistoManager.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(PrimaryGeneratorAction* kin)
:fPrimary(kin)
{
  fHistoManager = new HistoManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{ 
  delete fHistoManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{ 
  //initialize arrays
  //
  fDecayCount = fTimeCount = 0;
  for (G4int i=0; i<3; i++) fEkinTot[i] = fPbalance[i] = fEventTime[i] = 0. ;
  fPrimaryTime = 0.;
          
  //histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    analysisManager->OpenFile();
  }     
  
  //inform the runManager to save random number seed
  //
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);  
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

void RunAction::Balance(G4double Ekin, G4double Pbal)
{
  fDecayCount++;
  fEkinTot[0] += Ekin;
  //update min max  
  if (fDecayCount == 1) fEkinTot[1] = fEkinTot[2] = Ekin;
  if (Ekin < fEkinTot[1]) fEkinTot[1] = Ekin;
  if (Ekin > fEkinTot[2]) fEkinTot[2] = Ekin;
  
  fPbalance[0] += Pbal;
  //update min max   
  if (fDecayCount == 1) fPbalance[1] = fPbalance[2] = Pbal;  
  if (Pbal < fPbalance[1]) fPbalance[1] = Pbal;
  if (Pbal > fPbalance[2]) fPbalance[2] = Pbal;    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EventTiming(G4double time)
{
  fTimeCount++;  
  fEventTime[0] += time;
  if (fTimeCount == 1) fEventTime[1] = fEventTime[2] = time;  
  if (time < fEventTime[1]) fEventTime[1] = time;
  if (time > fEventTime[2]) fEventTime[2] = time;             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::PrimaryTiming(G4double ptime)
{
  fPrimaryTime += ptime;
}
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
 G4int nbEvents = run->GetNumberOfEvent();
 if (nbEvents == 0) { return; }
 
 G4ParticleDefinition* particle = fPrimary->GetParticleGun()
                                         ->GetParticleDefinition();
 G4String partName = particle->GetParticleName();
 G4double eprimary = fPrimary->GetParticleGun()->GetParticleEnergy();
 
 G4cout << "\n ======================== run summary ======================";  
 G4cout << "\n The run was " << nbEvents << " " << partName << " of "
        << G4BestUnit(eprimary,"Energy");
 G4cout << "\n ===========================================================\n";
 G4cout << G4endl;

 G4int prec = 4, wid = prec + 2;
 G4int dfprec = G4cout.precision(prec);
      
 //particle count
 //
 G4cout << " Nb of generated particles: \n" << G4endl;
     
 std::map<G4String,G4int>::iterator it;               
 for (it = fParticleCount.begin(); it != fParticleCount.end(); it++) { 
    G4String name = it->first;
    G4int count   = it->second;
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

 if (fDecayCount > 0) {
    G4double Ebmean = fEkinTot[0]/fDecayCount;
    G4double Pbmean = fPbalance[0]/fDecayCount;
         
    G4cout << "\n   Ekin Total (Q): mean = "
           << std::setw(wid) << G4BestUnit(Ebmean, "Energy")
           << "\t( "  << G4BestUnit(fEkinTot[1], "Energy")
           << " --> " << G4BestUnit(fEkinTot[2], "Energy")
           << ")" << G4endl;    
           
    G4cout << "\n   Momentum balance (excluding gamma desexcitation): mean = " 
           << std::setw(wid) << G4BestUnit(Pbmean, "Energy")
           << "\t( "  << G4BestUnit(fPbalance[1], "Energy")
           << " --> " << G4BestUnit(fPbalance[2], "Energy")
           << ")" << G4endl;
 }
            
 //total time of life
 //
 if (fTimeCount > 0) {
    G4double Tmean = fEventTime[0]/fTimeCount;
    G4double halfLife = Tmean*std::log(2.);
   
    G4cout << "\n   Total time of life : mean = "
           << std::setw(wid) << G4BestUnit(Tmean, "Time")
           << "  half-life = "
           << std::setw(wid) << G4BestUnit(halfLife, "Time")
           << "   ( "  << G4BestUnit(fEventTime[1], "Time")
           << " --> "  << G4BestUnit(fEventTime[2], "Time")
           << ")" << G4endl;
 }
            
 //activity of primary ion
 //
 G4double pTimeMean = fPrimaryTime/nbEvents;
 G4double molMass = particle->GetAtomicMass()*g/mole;
 G4double nAtoms_perUnitOfMass = Avogadro/molMass;
 G4double Activity_perUnitOfMass = 0.0;
 if (pTimeMean > 0.0)
   { Activity_perUnitOfMass = nAtoms_perUnitOfMass/pTimeMean; }
   
 G4cout << "\n   Activity of " << partName << " = "
            << std::setw(wid)  << Activity_perUnitOfMass*g/becquerel
            << " Bq/g   ("     << Activity_perUnitOfMass*g/curie
            << " Ci/g) \n" 
            << G4endl;
                                            
 // remove all contents in fParticleCount
 // 
 fParticleCount.clear(); 
 fEmean.clear();  fEmin.clear(); fEmax.clear();

 // restore default precision
 // 
 G4cout.precision(dfprec);
            
 //normalize and save histograms
 //
 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance(); 
 G4double factor = 100./nbEvents;
 analysisManager->ScaleH1(1,factor);
 analysisManager->ScaleH1(2,factor);
 analysisManager->ScaleH1(3,factor);
 analysisManager->ScaleH1(4,factor);
 analysisManager->ScaleH1(5,factor);   
 if ( analysisManager->IsActive() ) {
  analysisManager->Write();
  analysisManager->CloseFile();
 } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
