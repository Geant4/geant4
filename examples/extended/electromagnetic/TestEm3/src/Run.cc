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
/// \file electromagnetic/TestEm3/src/Run.cc
/// \brief Implementation of the Run class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"
#include "EmAcceptance.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det)
: G4Run(),
  fDetector(det), 
  fParticle(nullptr), fEkin(0.),
  fChargedStep(0), fNeutralStep(0),
  fN_gamma(0), fN_elec(0), fN_pos(0),
  fApplyLimit(false)
{
  //initialize cumulative quantities
  //
  for (G4int k=0; k<kMaxAbsor; k++) {
    fSumEAbs[k] = fSum2EAbs[k]  = fSumLAbs[k] = fSum2LAbs[k] = 0.;
    fEnergyDeposit[k].clear();
    fEdeptrue[k] = fRmstrue[k] = 1.;
    fLimittrue[k] = DBL_MAX;  
  }
  
  //initialize Eflow
  //
  G4int nbPlanes = (fDetector->GetNbOfLayers())*(fDetector->GetNbOfAbsor()) + 2;
  fEnergyFlow.resize(nbPlanes);
  fLateralEleak.resize(nbPlanes);
  for (G4int k=0; k<nbPlanes; k++) {fEnergyFlow[k] = fLateralEleak[k] = 0.; }  
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

void Run::FillPerEvent(G4int kAbs, G4double EAbs, G4double LAbs)
{
  //accumulate statistic with restriction
  //
  if(fApplyLimit) fEnergyDeposit[kAbs].push_back(EAbs);
  fSumEAbs[kAbs]  += EAbs;  fSum2EAbs[kAbs]  += EAbs*EAbs;
  fSumLAbs[kAbs]  += LAbs;  fSum2LAbs[kAbs]  += LAbs*LAbs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SumEnergyFlow(G4int plane, G4double Eflow)
{
  fEnergyFlow[plane] += Eflow;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
                                            
void Run::SumLateralEleak(G4int cell, G4double Eflow)
{
  fLateralEleak[cell] += Eflow;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
                                            
void Run::AddChargedStep() 
{
  fChargedStep += 1.0; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
void Run::AddNeutralStep() 
{
  fNeutralStep += 1.0; 
}
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddSecondaryTrack(const G4Track* track)
{
  const G4ParticleDefinition* d = track->GetDefinition();
  if(d == G4Gamma::Gamma()) { ++fN_gamma; }
  else if (d == G4Electron::Electron()) { ++fN_elec; }
  else if (d == G4Positron::Positron()) { ++fN_pos; }
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
  for (G4int k=0; k<kMaxAbsor; k++) {
    fSumEAbs[k]  += localRun->fSumEAbs[k]; 
    fSum2EAbs[k] += localRun->fSum2EAbs[k]; 
    fSumLAbs[k]  += localRun->fSumLAbs[k]; 
    fSum2LAbs[k] += localRun->fSum2LAbs[k];
  }
   
  G4int nbPlanes = (fDetector->GetNbOfLayers())*(fDetector->GetNbOfAbsor()) + 2;
  for (G4int k=0; k<nbPlanes; k++) {
    fEnergyFlow[k]   += localRun->fEnergyFlow[k];
    fLateralEleak[k] += localRun->fLateralEleak[k];
  }
  
       
  fChargedStep += localRun->fChargedStep;  
  fNeutralStep += localRun->fNeutralStep;
  
  fN_gamma += localRun->fN_gamma;    
  fN_elec  += localRun->fN_elec;
  fN_pos   += localRun->fN_pos;
     
  fApplyLimit = localRun->fApplyLimit;
  
  for (G4int k=0; k<kMaxAbsor; k++) {
    fEdeptrue[k]  = localRun->fEdeptrue[k]; 
    fRmstrue[k]   = localRun->fRmstrue[k]; 
    fLimittrue[k] = localRun->fLimittrue[k]; 
  }
    
  G4Run::Merge(run); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun()
{
  G4int nEvt = numberOfEvent;
  G4double  norm = G4double(nEvt);
  if(norm > 0) norm = 1./norm;
  G4double qnorm = std::sqrt(norm);

  fChargedStep *= norm;
  fNeutralStep *= norm;

  //compute and print statistic
  //
  G4double beamEnergy = fEkin;
  G4double sqbeam = std::sqrt(beamEnergy/GeV);

  G4double MeanEAbs,MeanEAbs2,rmsEAbs,resolution,rmsres;
  G4double MeanLAbs,MeanLAbs2,rmsLAbs;

  std::ios::fmtflags mode = G4cout.flags();
  G4int  prec = G4cout.precision(2);
  G4cout << "\n------------------------------------------------------------\n";
  G4cout << std::setw(14) << "material"
         << std::setw(17) << "Edep       RMS"
         << std::setw(33) << "sqrt(E0(GeV))*rmsE/Emean"
         << std::setw(23) << "total tracklen \n \n";

  for (G4int k=1; k<=fDetector->GetNbOfAbsor(); k++)
    {
      MeanEAbs  = fSumEAbs[k]*norm;
      MeanEAbs2 = fSum2EAbs[k]*norm;
      rmsEAbs  = std::sqrt(std::abs(MeanEAbs2 - MeanEAbs*MeanEAbs));
      //G4cout << "k= " << k << "  RMS= " <<  rmsEAbs 
      //     << "  fApplyLimit: " << fApplyLimit << G4endl;
      if(fApplyLimit) {
        G4int    nn    = 0;
        G4double sume  = 0.0;
        G4double sume2 = 0.0;
        // compute trancated means  
        G4double lim   = rmsEAbs * 2.5;
        for(G4int i=0; i<nEvt; i++) {
          G4double e = (fEnergyDeposit[k])[i];
          if(std::abs(e - MeanEAbs) < lim) {
            sume  += e;
            sume2 += e*e;
            nn++;
          }
        }
        G4double norm1 = G4double(nn);
        if(norm1 > 0.0) norm1 = 1.0/norm1;
        MeanEAbs  = sume*norm1;
        MeanEAbs2 = sume2*norm1;
        rmsEAbs  = std::sqrt(std::abs(MeanEAbs2 - MeanEAbs*MeanEAbs));
      }

      resolution= 100.*sqbeam*rmsEAbs/MeanEAbs;
      rmsres    = resolution*qnorm;

      // Save mean and RMS
      fSumEAbs[k] = MeanEAbs;
      fSum2EAbs[k] = rmsEAbs;

      MeanLAbs  = fSumLAbs[k]*norm;
      MeanLAbs2 = fSum2LAbs[k]*norm;
      rmsLAbs  = std::sqrt(std::abs(MeanLAbs2 - MeanLAbs*MeanLAbs));

      //print
      //
      G4cout
       << std::setw(14) << fDetector->GetAbsorMaterial(k)->GetName() << ": "
       << std::setprecision(5)
       << std::setw(6) << G4BestUnit(MeanEAbs,"Energy") << " :  "
       << std::setprecision(4)
       << std::setw(5) << G4BestUnit( rmsEAbs,"Energy")  
       << std::setw(10) << resolution  << " +- " 
       << std::setw(5) << rmsres << " %"
       << std::setprecision(3)
       << std::setw(10) << G4BestUnit(MeanLAbs,"Length")  << " +- "
       << std::setw(4) << G4BestUnit( rmsLAbs,"Length")
       << G4endl;
    }
  G4cout << "\n------------------------------------------------------------\n";

  G4cout << " Beam particle " 
         << fParticle->GetParticleName()
         << "  E = " << G4BestUnit(beamEnergy,"Energy") << G4endl;
  G4cout << " Mean number of gamma       " << (G4double)fN_gamma*norm << G4endl;
  G4cout << " Mean number of e-          " << (G4double)fN_elec*norm << G4endl;
  G4cout << " Mean number of e+          " << (G4double)fN_pos*norm << G4endl;
  G4cout << std::setprecision(6)
         << " Mean number of charged steps  " << fChargedStep << G4endl;
  G4cout << " Mean number of neutral steps  " << fNeutralStep << G4endl;
  G4cout << "------------------------------------------------------------\n";
  
  //Energy flow
  //
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  G4int Idmax = (fDetector->GetNbOfLayers())*(fDetector->GetNbOfAbsor());
  for (G4int Id=1; Id<=Idmax+1; Id++) {
    analysis->FillH1(2*kMaxAbsor+1, (G4double)Id, fEnergyFlow[Id]);
    analysis->FillH1(2*kMaxAbsor+2, (G4double)Id, fLateralEleak[Id]);
  }
  
  //Energy deposit from energy flow balance
  //
  G4double EdepTot[kMaxAbsor];
  for (G4int k=0; k<kMaxAbsor; k++) EdepTot[k] = 0.;
  
  G4int nbOfAbsor = fDetector->GetNbOfAbsor();
  for (G4int Id=1; Id<=Idmax; Id++) {
    G4int iAbsor = Id%nbOfAbsor; if (iAbsor==0) iAbsor = nbOfAbsor;
    EdepTot[iAbsor] += (fEnergyFlow[Id]-fEnergyFlow[Id+1]-fLateralEleak[Id]);
  }
  
  G4cout << std::setprecision(3)
         << "\n Energy deposition from Energy flow balance : \n"
         << std::setw(10) << "  material \t Total Edep \n \n";
  G4cout.precision(6);
  
  for (G4int k=1; k<=nbOfAbsor; k++) {
    EdepTot [k] *= norm;
    G4cout << std::setw(10) << fDetector->GetAbsorMaterial(k)->GetName() << ":"
           << "\t " << G4BestUnit(EdepTot [k],"Energy") << "\n";
  }
  
  G4cout << "\n------------------------------------------------------------\n" 
         << G4endl;

  // Acceptance
  EmAcceptance acc;
  G4bool isStarted = false;
  for (G4int j=1; j<=fDetector->GetNbOfAbsor(); j++) {
    if (fLimittrue[j] < DBL_MAX) {
      if (!isStarted) {
        acc.BeginOfAcceptance("Sampling Calorimeter",nEvt);
        isStarted = true;
      }
      MeanEAbs = fSumEAbs[j];
      rmsEAbs  = fSum2EAbs[j];
      G4String mat = fDetector->GetAbsorMaterial(j)->GetName();
      acc.EmAcceptanceGauss("Edep"+mat, nEvt, MeanEAbs,
                             fEdeptrue[j], fRmstrue[j], fLimittrue[j]);
      acc.EmAcceptanceGauss("Erms"+mat, nEvt, rmsEAbs,
                             fRmstrue[j], fRmstrue[j], 2.0*fLimittrue[j]);
    }
  }
  if(isStarted) acc.EndOfAcceptance();

  //normalize histograms
  //
  for (G4int ih = kMaxAbsor+1; ih < kMaxHisto; ih++) {
    analysis->ScaleH1(ih,norm/MeV);
  }
  
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetEdepAndRMS(G4int i, G4double edep, G4double rms, G4double lim)
{
  if (i>=0 && i<kMaxAbsor) {
    fEdeptrue [i] = edep;
    fRmstrue  [i] = rms;
    fLimittrue[i] = lim;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetApplyLimit(G4bool val)
{
  fApplyLimit = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
