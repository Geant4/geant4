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
/// \file electromagnetic/TestEm11/src/Run.cc
/// \brief Implementation of the Run class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "HistoManager.hh"

#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det)
: G4Run(),
  fDetector(det), 
  fParticle(0), fEkin(0.),
  nbOfModules(0), nbOfLayers(0), kLayerMax(0),
  EtotCalor(0.), Etot2Calor(0.), EvisCalor(0.), Evis2Calor(0.),
  Eleak(0.), Eleak2(0.)
{
  nbOfModules = fDetector->GetNbModules();	 	
  nbOfLayers  = fDetector->GetNbLayers();
  kLayerMax = nbOfModules*nbOfLayers + 1;
  
  //initialize vectors
  //
  EtotLayer.resize(kLayerMax); Etot2Layer.resize(kLayerMax);
  EvisLayer.resize(kLayerMax); Evis2Layer.resize(kLayerMax);			
  for (G4int k=0; k<kLayerMax; k++) {
    EtotLayer[k] = Etot2Layer[k] = EvisLayer[k] = Evis2Layer[k] = 0.0;
  }
  
  EtotCalor = Etot2Calor = EvisCalor = Evis2Calor = Eleak = Eleak2 = 0.;
  EdLeak[0] = EdLeak[1] = EdLeak[2] = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SumEvents_1(G4int layer, G4double Etot, G4double Evis)
{
  //accumulate statistic per layer
  //
  EtotLayer[layer] += Etot;  Etot2Layer[layer] += Etot*Etot;
  EvisLayer[layer] += Evis;  Evis2Layer[layer] += Evis*Evis;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SumEvents_2(G4double etot, G4double evis, G4double eleak)
{
  //accumulate statistic for full calorimeter
  //
  EtotCalor += etot;  Etot2Calor += etot*etot; 	
  EvisCalor += evis;  Evis2Calor += evis*evis; 
  Eleak += eleak;  Eleak2 += eleak*eleak;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::DetailedLeakage(G4int icase, G4double energy)
{
  //forward, backward, lateral leakage
  //
  EdLeak[icase] += energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);

  // pass information about primary particle
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;

  // accumulate sums
  //
  for (G4int k=0; k<kLayerMax; k++) {
    EtotLayer[k]  += localRun->EtotLayer[k];
    Etot2Layer[k] += localRun->Etot2Layer[k];	 
    EvisLayer[k]  += localRun->EvisLayer[k];
    Evis2Layer[k] += localRun->Evis2Layer[k];	 
  }
  
  EtotCalor  += localRun->EtotCalor;  
  Etot2Calor += localRun->Etot2Calor;  
  EvisCalor  += localRun->EvisCalor;
  Evis2Calor += localRun->Evis2Calor;   
  Eleak      += localRun->Eleak;  
  Eleak2     += localRun->Eleak2;
  EdLeak[0]  += localRun->EdLeak[0];
  EdLeak[1]  += localRun->EdLeak[1];
  EdLeak[2]  += localRun->EdLeak[2];
  
  G4Run::Merge(run); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun()
{
  //calorimeter
  //
  fDetector->PrintCalorParameters();
 
  //run conditions
  //   
  G4String partName = fParticle->GetParticleName();
  G4int nbEvents = numberOfEvent;

  G4int prec = G4cout.precision(3);

  G4cout << " The run was " << nbEvents << " " << partName << " of "
         << G4BestUnit(fEkin,"Energy") << " through the calorimeter" << G4endl;
 
  G4cout << "------------------------------------------------------------"
         << G4endl;
   
  //if no events, return
  //
  if (nbEvents == 0) return;

  //compute and print statistic
  //
  std::ios::fmtflags mode = G4cout.flags(); 
   
  // energy in layers
  //
  G4cout.precision(prec);	 
  G4cout << "\n             " 
         << "total Energy          (rms/mean)      "
         << "visible Energy        (rms/mean)" << G4endl;
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  G4double meanEtot,meanEtot2,varianceEtot,rmsEtot,resEtot;  
  G4double meanEvis,meanEvis2,varianceEvis,rmsEvis,resEvis;
  
  for (G4int i1=1; i1<kLayerMax; i1++) {
    //total energy
    meanEtot  = EtotLayer[i1] /nbEvents;
    meanEtot2 = Etot2Layer[i1]/nbEvents;    
    varianceEtot = meanEtot2 - meanEtot*meanEtot;
    resEtot = rmsEtot = 0.;
    if (varianceEtot > 0.) rmsEtot = std::sqrt(varianceEtot);
    if (meanEtot > 0.) resEtot = 100*rmsEtot/meanEtot;
    analysisManager->FillH1(3, i1+0.5, meanEtot);
   	  
    //visible energy
    meanEvis  = EvisLayer[i1] /nbEvents;
    meanEvis2 = Evis2Layer[i1]/nbEvents;    
    varianceEvis = meanEvis2 - meanEvis*meanEvis;
    resEvis = rmsEvis = 0.;
    if (varianceEvis > 0.) rmsEvis = std::sqrt(varianceEvis);
  	if (meanEvis > 0.) resEvis = 100*rmsEvis/meanEvis;
    analysisManager->FillH1(4, i1+0.5, meanEvis);

    //print
    //
    G4cout
      << "\n   layer " << i1 << ": "
      << std::setprecision(5)
      << std::setw(6) << G4BestUnit(meanEtot,"Energy") << " +- "
      << std::setprecision(4)
      << std::setw(5) << G4BestUnit( rmsEtot,"Energy") << "  ("
      << std::setprecision(2) 
      << std::setw(3) << resEtot  << " %)" 
      << "     "
      << std::setprecision(5)
      << std::setw(6) << G4BestUnit(meanEvis,"Energy") << " +- "
      << std::setprecision(4)
      << std::setw(5) << G4BestUnit( rmsEvis,"Energy") << "  ("
      << std::setprecision(2) 
      << std::setw(3) << resEvis  << " %)"; 
  }
  G4cout << G4endl;

  //calorimeter: total energy
  meanEtot  = EtotCalor /nbEvents;
  meanEtot2 = Etot2Calor/nbEvents;
  varianceEtot = meanEtot2 - meanEtot*meanEtot;
  resEtot = rmsEtot = 0.;
  if (varianceEtot > 0.) rmsEtot = std::sqrt(varianceEtot);
  if (meanEtot > 0.) resEtot = 100*rmsEtot/meanEtot;

  //calorimeter: visible energy
  meanEvis  = EvisCalor /nbEvents;
  meanEvis2 = Evis2Calor/nbEvents;
  varianceEvis = meanEvis2 - meanEvis*meanEvis;
  resEvis = rmsEvis = 0.;
  if (varianceEvis > 0.) rmsEvis = std::sqrt(varianceEvis);
  if (meanEvis > 0.) resEvis = 100*rmsEvis/meanEvis;
      
  //print
  //
  G4cout
    << "\n   total calor : "
    << std::setprecision(5)
    << std::setw(6) << G4BestUnit(meanEtot,"Energy") << " +- "
    << std::setprecision(4)
    << std::setw(5) << G4BestUnit( rmsEtot,"Energy") << "  ("
    << std::setprecision(2) 
    << std::setw(3) << resEtot  << " %)" 
    << "     "
    << std::setprecision(5)
    << std::setw(6) << G4BestUnit(meanEvis,"Energy") << " +- "
    << std::setprecision(4)
    << std::setw(5) << G4BestUnit( rmsEvis,"Energy") << "  ("
    << std::setprecision(2) 
    << std::setw(3) << resEvis  << " %)";
                     
  G4cout << "\n------------------------------------------------------------"
         << G4endl;

  //leakage
  G4double meanEleak,meanEleak2,varianceEleak,rmsEleak,ratio;
  meanEleak  = Eleak /nbEvents;
  meanEleak2 = Eleak2/nbEvents;
  varianceEleak = meanEleak2 - meanEleak*meanEleak;
  rmsEleak = 0.;
  if (varianceEleak > 0.) rmsEleak = std::sqrt(varianceEleak);
  ratio = 100*meanEleak/fEkin;

  G4double forward = 100*EdLeak[0]/(nbEvents*fEkin);
  G4double bakward = 100*EdLeak[1]/(nbEvents*fEkin);
  G4double lateral = 100*EdLeak[2]/(nbEvents*fEkin);
      
  //print
  //
  G4cout
    << "\n   Leakage : "
    << std::setprecision(5)
    << std::setw(6) << G4BestUnit(meanEleak,"Energy") << " +- "
    << std::setprecision(4)
    << std::setw(5) << G4BestUnit( rmsEleak,"Energy") 
    << "\n   Eleak/Ebeam ="
    << std::setprecision(3) 
    << std::setw(4) << ratio  << " %  ( forward ="
    << std::setw(4) << forward  << " %;   backward ="
    << std::setw(4) << bakward  << " %;   lateral ="
    << std::setw(4) << lateral  << " %)"             
    << G4endl;
  
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);
  
  //normalize histograms
  G4double factor = 1./nbEvents;
  analysisManager->ScaleH1(5,factor);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
