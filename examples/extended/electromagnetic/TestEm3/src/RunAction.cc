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
/// \file electromagnetic/TestEm3/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "PrimaryGeneratorAction.hh"
#include "RunActionMessenger.hh"
#include "HistoManager.hh"
#include "EmAcceptance.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ProductionCutsTable.hh"
#include "G4LossTableManager.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* prim)
:fDetector(det), fPrimary(prim), fRunMessenger(0), fHistoManager(0)
{
  fRunMessenger = new RunActionMessenger(this);
  fHistoManager = new HistoManager();
  fApplyLimit = false;

  fChargedStep = fNeutralStep = 0.0;

  for (G4int k=0; k<MaxAbsor; k++) { fEdeptrue[k] = fRmstrue[k] = 1.;
                                    fLimittrue[k] = DBL_MAX;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete fRunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  // save Rndm status
  //
  ////G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  CLHEP::HepRandom::showEngineStatus();

  //initialize cumulative quantities
  //
  for (G4int k=0; k<MaxAbsor; k++) {
    fSumEAbs[k] = fSum2EAbs[k]  = fSumLAbs[k] = fSum2LAbs[k] = 0.;
    fEnergyDeposit[k].clear();  
  }

  fChargedStep = fNeutralStep = 0.0;

  fN_gamma = 0;
  fN_elec  = 0;
  fN_pos   = 0;

  //initialize Eflow
  //
  G4int nbPlanes = (fDetector->GetNbOfLayers())*(fDetector->GetNbOfAbsor()) + 2;
  fEnergyFlow.resize(nbPlanes);
  fLateralEleak.resize(nbPlanes);
  for (G4int k=0; k<nbPlanes; k++) {fEnergyFlow[k] = fLateralEleak[k] = 0.; }
  
  //histograms
  //
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  if (analysis->IsActive()) analysis->OpenFile();
   
  //example of print dEdx tables
  //
  ////PrintDedxTables();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::FillPerEvent(G4int kAbs, G4double EAbs, G4double LAbs)
{
  //accumulate statistic with restriction
  //
  if(fApplyLimit) fEnergyDeposit[kAbs].push_back(EAbs);
  fSumEAbs[kAbs]  += EAbs;  fSum2EAbs[kAbs]  += EAbs*EAbs;
  fSumLAbs[kAbs]  += LAbs;  fSum2LAbs[kAbs]  += LAbs*LAbs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int nEvt = aRun->GetNumberOfEvent();
  G4double  norm = G4double(nEvt);
  if(norm > 0) norm = 1./norm;
  G4double qnorm = std::sqrt(norm);

  fChargedStep *= norm;
  fNeutralStep *= norm;

  //compute and print statistic
  //
  G4double beamEnergy = fPrimary->GetParticleGun()->GetParticleEnergy();
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
         << fPrimary->GetParticleGun()->
    GetParticleDefinition()->GetParticleName()
         << "  E = " << G4BestUnit(beamEnergy,"Energy") << G4endl;
  G4cout << " Mean number of gamma          " << (G4double)fN_gamma*norm << G4endl;
  G4cout << " Mean number of e-             " << (G4double)fN_elec*norm << G4endl;
  G4cout << " Mean number of e+             " << (G4double)fN_pos*norm << G4endl;
  G4cout << std::setprecision(6)
	 << " Mean number of charged steps  " << fChargedStep << G4endl;
  G4cout << " Mean number of neutral steps  " << fNeutralStep << G4endl;
  G4cout << "------------------------------------------------------------\n";
  
  //Energy flow
  //
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  G4int Idmax = (fDetector->GetNbOfLayers())*(fDetector->GetNbOfAbsor());
  for (G4int Id=1; Id<=Idmax+1; Id++) {
    analysis->FillH1(2*MaxAbsor+1, (G4double)Id, fEnergyFlow[Id]);
    analysis->FillH1(2*MaxAbsor+2, (G4double)Id, fLateralEleak[Id]);
  }
  
  //Energy deposit from energy flow balance
  //
  G4double EdepTot[MaxAbsor];
  for (G4int k=0; k<MaxAbsor; k++) EdepTot[k] = 0.;
  
  G4int nbOfAbsor = fDetector->GetNbOfAbsor();
  for (G4int Id=1; Id<=Idmax; Id++) {
    G4int iAbsor = Id%nbOfAbsor; if (iAbsor==0) iAbsor = nbOfAbsor;
    EdepTot [iAbsor] += (fEnergyFlow[Id] - fEnergyFlow[Id+1] - fLateralEleak[Id]);
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
    
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);

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
  for (G4int ih = MaxAbsor+1; ih < MaxHisto; ih++) {
    analysis->ScaleH1(ih,norm/MeV);
  }
  
  //save histograms
  if (analysis->IsActive()) {   
    analysis->Write();
    analysis->CloseFile();
  }    

  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::PrintDedxTables()
{
  //Print dE/dx tables with binning identical to the Geant3 JMATE bank.
  //The printout is readable as Geant3 ffread data cards (by the program g4mat).
  //
  const G4double tkmin=10*keV, tkmax=10*TeV;
  const G4int nbin=90;
  G4double tk[nbin];

  const G4int ncolumn = 5;

  //compute the kinetic energies
  //
  const G4double dp = std::log10(tkmax/tkmin)/nbin;
  const G4double dt = std::pow(10.,dp);
  tk[0] = tkmin;
  for (G4int i=1; i<nbin; ++i) tk[i] = tk[i-1]*dt;

  //print the kinetic energies
  //
  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  G4int  prec = G4cout.precision(3);

  G4cout << "\n kinetic energies \n ";
  for (G4int j=0; j<nbin; ++j) {
    G4cout << G4BestUnit(tk[j],"Energy") << "\t";
    if ((j+1)%ncolumn == 0) G4cout << "\n ";
  }
  G4cout << G4endl;

  //print the dE/dx tables
  //
  G4cout.setf(std::ios::scientific,std::ios::floatfield);

  G4ParticleDefinition*
  part = fPrimary->GetParticleGun()->GetParticleDefinition();
  
  G4ProductionCutsTable* theCoupleTable =
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();
  const G4MaterialCutsCouple* couple = 0;

  for (G4int iab=1;iab <= fDetector->GetNbOfAbsor(); iab++)
     {
      G4Material* mat = fDetector->GetAbsorMaterial(iab);
      G4int index = 0;
      for (size_t i=0; i<numOfCouples; i++) {
         couple = theCoupleTable->GetMaterialCutsCouple(i);
         if (couple->GetMaterial() == mat) {index = i; break;}
      }
      G4cout << "\nLIST";
      G4cout << "\nC \nC  dE/dx (MeV/cm) for " << part->GetParticleName()
             << " in " << mat ->GetName() << "\nC";
      G4cout << "\nKINE   (" << part->GetParticleName() << ")";
      G4cout << "\nMATE   (" << mat ->GetName() << ")";
      G4cout.precision(2);
      G4cout << "\nERAN  " << tkmin/GeV << " (ekmin)\t"
                           << tkmax/GeV << " (ekmax)\t"
                           << nbin      << " (nekbin)";
      G4double cutgam =
         (*(theCoupleTable->GetEnergyCutsVector(idxG4GammaCut)))[index];
      if (cutgam < tkmin) cutgam = tkmin;
      if (cutgam > tkmax) cutgam = tkmax;
      G4double cutele =
         (*(theCoupleTable->GetEnergyCutsVector(idxG4ElectronCut)))[index];
      if (cutele < tkmin) cutele = tkmin;
      if (cutele > tkmax) cutele = tkmax;
      G4cout << "\nCUTS  " << cutgam/GeV << " (cutgam)\t"
                           << cutele/GeV << " (cutele)";

      G4cout.precision(6);
      G4cout << "\nG4VAL \n ";
      for (G4int l=0;l<nbin; ++l)
         {
           G4double dedx = G4LossTableManager::Instance()
                                               ->GetDEDX(part,tk[l],couple);
           G4cout << dedx/(MeV/cm) << "\t";
           if ((l+1)%ncolumn == 0) G4cout << "\n ";
         }
      G4cout << G4endl;
     }

  G4cout.precision(prec);
  G4cout.setf(mode,std::ios::floatfield);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddSecondaryTrack(const G4Track* track)
{
  const G4ParticleDefinition* d = track->GetDefinition();
  if(d == G4Gamma::Gamma()) { ++fN_gamma; }
  else if (d == G4Electron::Electron()) { ++fN_elec; }
  else if (d == G4Positron::Positron()) { ++fN_pos; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetEdepAndRMS(G4int i, G4double edep, G4double rms, G4double lim)
{
  if (i>=0 && i<MaxAbsor) {
    fEdeptrue [i] = edep;
    fRmstrue  [i] = rms;
    fLimittrue[i] = lim;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
