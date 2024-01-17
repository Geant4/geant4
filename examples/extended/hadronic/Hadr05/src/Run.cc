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
/// \file Run.cc
/// \brief Implementation of the Run class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4ParticleDefinition.hh"
#include "G4Track.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det)
: fDetector(det)
{
  //initialize energy deposited per absorber
  //
  for (G4int k=0; k<kMaxAbsor; k++) {
    fSumEAbs[k] = fSum2EAbs[k]  = fSumLAbs[k] = fSum2LAbs[k] = 0.;
  }
  
  // initialize total energy deposited
  //
  fEdepTot = fEdepTot2 = 0.;
  
  // initialize leakage
  //
  fEnergyLeak[0] = fEnergyLeak[1] = 0.;
  fEleakTot = fEleakTot2 = 0.;
  
  // initialize total energy released
  //
  fEtotal = fEtotal2 = 0.;
      
  //initialize Eflow
  //
  G4int nbPlanes = (fDetector->GetNbOfLayers())*(fDetector->GetNbOfAbsor()) + 2;
  fEnergyFlow.resize(nbPlanes);
  for (G4int k=0; k<nbPlanes; k++) {fEnergyFlow[k] = 0.; }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin = energy;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::CountProcesses(const G4VProcess* process) 
{
  if (process == nullptr) return;
  G4String procName = process->GetProcessName();
  std::map<G4String,G4int>::iterator it = fProcCounter.find(procName);
  if ( it == fProcCounter.end()) {
    fProcCounter[procName] = 1;
  }
  else {
    fProcCounter[procName]++; 
  }
}
   
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SumEdepPerAbsorber(G4int kAbs, G4double EAbs, G4double LAbs)
{
  //accumulate statistic with restriction
  //
  fSumEAbs[kAbs]  += EAbs;  fSum2EAbs[kAbs]  += EAbs*EAbs;
  fSumLAbs[kAbs]  += LAbs;  fSum2LAbs[kAbs]  += LAbs*LAbs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SumEnergies(G4double edeptot, G4double eleak0, G4double eleak1)
{
  fEdepTot += edeptot; fEdepTot2 += edeptot*edeptot;
 
  fEnergyLeak[0] += eleak0; fEnergyLeak[1] += eleak1;
  G4double eleaktot = eleak0 + eleak1;
  fEleakTot += eleaktot; fEleakTot2 += eleaktot*eleaktot;
  
  G4double etotal = edeptot + eleaktot;
  fEtotal += etotal; fEtotal2 += etotal*etotal;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SumEnergyFlow(G4int plane, G4double Eflow)
{
  fEnergyFlow[plane] += Eflow;
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
  
  fEdepTot  += localRun->fEdepTot;
  fEdepTot2 += localRun->fEdepTot2;
  
  fEnergyLeak[0]  += localRun->fEnergyLeak[0];
  fEnergyLeak[1]  += localRun->fEnergyLeak[1];
  
  fEleakTot  += localRun->fEleakTot;
  fEleakTot2 += localRun->fEleakTot2;
  
  fEtotal  += localRun->fEtotal;
  fEtotal2 += localRun->fEtotal2;  
     
  G4int nbPlanes = (fDetector->GetNbOfLayers())*(fDetector->GetNbOfAbsor()) + 2;
  for (G4int k=0; k<nbPlanes; k++) {
    fEnergyFlow[k]   += localRun->fEnergyFlow[k];
  }
  
  //map: processes count
  std::map<G4String,G4int>::const_iterator itp;
  for ( itp = localRun->fProcCounter.begin();
        itp != localRun->fProcCounter.end(); ++itp ) {

    G4String procName = itp->first;
    G4int localCount = itp->second;
    if ( fProcCounter.find(procName) == fProcCounter.end()) {
      fProcCounter[procName] = localCount;
    }
    else {
      fProcCounter[procName] += localCount;
    }  
  }
      
  G4Run::Merge(run); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun()
{
  //run condition
  //   
  G4String Particle = fParticle->GetParticleName();    
  G4cout << "\n ---> The run is " << numberOfEvent << " "<< Particle << " of "
         << G4BestUnit(fEkin,"Energy") << " through calorimeter" << G4endl;
         
  //frequency of processes
  //
  G4cout << "\n Process calls frequency :" << G4endl;
  G4int index = 0;
  std::map<G4String,G4int>::iterator it;    
  for (it = fProcCounter.begin(); it != fProcCounter.end(); it++) {
     G4String procName = it->first;
     G4int    count    = it->second;
     G4String space = " "; if (++index%3 == 0) space = "\n";
     G4cout << " " << std::setw(22) << procName << "="<< std::setw(10) << count
            << space;
  }
  
  G4cout << G4endl;
  G4int nEvt = numberOfEvent;
  G4double  norm = G4double(nEvt);
  if(norm > 0) norm = 1./norm;
  G4double qnorm = std::sqrt(norm);

  //energy deposit per absorber
  //
  G4double beamEnergy = fEkin;
  G4double sqbeam = std::sqrt(beamEnergy/GeV);

  G4double MeanEAbs,MeanEAbs2,rmsEAbs,resolution,rmsres;
  G4double MeanLAbs,MeanLAbs2,rmsLAbs;

  std::ios::fmtflags mode = G4cout.flags();
  G4int  prec = G4cout.precision(2);
  G4cout << "\n------------------------------------------------------------\n";
  G4cout << std::setw(16) << "material"
         << std::setw(22) << "Edep        rmsE"
         << std::setw(31) << "sqrt(E0(GeV))*rmsE/Edep"
         << std::setw(23) << "total tracklen \n \n";

  for (G4int k=1; k<=fDetector->GetNbOfAbsor(); k++)
    {
      MeanEAbs  = fSumEAbs[k]*norm;
      MeanEAbs2 = fSum2EAbs[k]*norm;
      rmsEAbs  = std::sqrt(std::abs(MeanEAbs2 - MeanEAbs*MeanEAbs));

      resolution= 100.*sqbeam*rmsEAbs/MeanEAbs;
      rmsres    = resolution*qnorm;

      MeanLAbs  = fSumLAbs[k]*norm;
      MeanLAbs2 = fSum2LAbs[k]*norm;
      rmsLAbs  = std::sqrt(std::abs(MeanLAbs2 - MeanLAbs*MeanLAbs));

      //print
      //
      G4cout
       << std::setw(2) << k
       << std::setw(14) << fDetector->GetAbsorMaterial(k)->GetName()
       << std::setprecision(5)
       << std::setw(10) << G4BestUnit(MeanEAbs,"Energy")
       << std::setprecision(4)
       << std::setw(8) << G4BestUnit( rmsEAbs,"Energy")  
       << std::setw(10) << resolution  << " +- "
       << std::setprecision(3)
       << std::setw(5) << rmsres << " %"
       << std::setprecision(4)
       << std::setw(12) << G4BestUnit(MeanLAbs,"Length")  << " +- "
       << std::setprecision(3)
       << std::setw(5) << G4BestUnit( rmsLAbs,"Length")
       << G4endl;
    }
 
  //total energy deposited
  //
  fEdepTot      /= nEvt;
  fEdepTot2     /= nEvt;
  G4double rmsEdep = std::sqrt(std::abs(fEdepTot2 - fEdepTot*fEdepTot));
  
  G4cout << "\n Total energy deposited = " << std::setprecision(4)
         << G4BestUnit(fEdepTot,"Energy")
	 << " +- " << G4BestUnit(rmsEdep, "Energy") << G4endl;	   
         
  //Energy leakage
  //
  fEnergyLeak[0] /= nEvt;  
  fEnergyLeak[1] /= nEvt;
  fEleakTot      /= nEvt;
  fEleakTot2     /= nEvt;
  G4double rmsEleak = std::sqrt(std::abs(fEleakTot2 - fEleakTot*fEleakTot));
  
  G4cout << " Leakage :  primary = "
         << G4BestUnit(fEnergyLeak[0],"Energy")
         << "   secondaries = "
         << G4BestUnit(fEnergyLeak[1],"Energy")
         << "  ---> total = " << G4BestUnit(fEleakTot, "Energy")
	 << " +- " << G4BestUnit(rmsEleak, "Energy") << G4endl;
	 
  //total energy released
  //
  fEtotal      /= nEvt;
  fEtotal2     /= nEvt;
  G4double rmsEtotal = std::sqrt(std::abs(fEtotal2 - fEtotal*fEtotal));
         
  G4cout << " Total energy released :  Edep + Eleak = "
         << G4BestUnit(fEtotal,"Energy")
	 << " +- " << G4BestUnit(rmsEtotal, "Energy") << G4endl;	            
  G4cout << "------------------------------------------------------------\n";
                     
  //Energy flow
  //
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  G4int Idmax = (fDetector->GetNbOfLayers())*(fDetector->GetNbOfAbsor());
  for (G4int Id=1; Id<=Idmax+1; Id++) {
    analysis->FillH1(2*kMaxAbsor+1, (G4double)Id, fEnergyFlow[Id]/nEvt);
  }
  
  //normalize histograms
  //
  for (G4int ih = kMaxAbsor+1; ih < 2*kMaxAbsor+1; ih++) {
    analysis->ScaleH1(ih,norm/MeV);
  }
  
  //remove all contents in fProcCounter
  fProcCounter.clear();
    
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
