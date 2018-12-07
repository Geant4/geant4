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
/// \file electromagnetic/TestEm8/src/Run.cc
/// \brief Implementation of the Run class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "G4Step.hh"
#include "G4Run.hh"
#include "G4LossTableManager.hh"
#include "G4ElectronIonPair.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "TestParameters.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run()
 : G4Run(), fElIonPair(0), fParam(TestParameters::GetPointer())
{
  fElIonPair = G4LossTableManager::Instance()->ElectronIonPair();
  fTotStepGas = fTotCluster = fMeanCluster = fOverflow = fTotEdep 
    = fStepGas = fCluster = fMaxEnergy = 0.0; 
  fEvt = fNbins = 0;
  fFactorALICE = fParam->GetFactorALICE();
  fWidthALICE = fParam->GetEnergySmear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::BeginOfRun()
{
  // initilise scoring
  fTotStepGas = fTotCluster = fMeanCluster = fOverflow = fTotEdep 
    = fStepGas = fCluster = 0.0; 
  fEvt = 0;

  fFactorALICE = fParam->GetFactorALICE();
  fWidthALICE = fParam->GetEnergySmear();
 
  SetVerbose(1);

  fNbins = fParam->GetNumberBins();
  fMaxEnergy = fParam->GetMaxEnergy();
  
  fEgas.resize(fNbins,0.0);
  fEdep.reset();

  if(fVerbose > 0) {
    G4int binsCluster = fParam->GetNumberBinsCluster();
    G4cout << " BinsCluster= " << binsCluster << "    BinsE= " <<  fNbins
           << "   Emax(keV)= " << fMaxEnergy/keV << G4endl;
    G4cout << " WidthALICE(keV)= " << fWidthALICE/keV 
           << "      FactorALICE= " << fFactorALICE << G4endl;

  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun()
{
  G4int nEvt = GetNumberOfEvent();
  G4double norm = (nEvt > 0) ? 1.0/(G4double)nEvt : 0.0; 

  fTotStepGas  *= norm;
  fTotCluster  *= norm;
  fMeanCluster *= norm;
  fOverflow    *= norm;

  G4double y1 = fEdep.mean();
  G4double y2 = fEdep.rms();

  G4double de = fMaxEnergy/G4double(fNbins);  
  G4double x1 = -de*0.5; 

  fFactorALICE = fParam->GetFactorALICE();

  G4cout << " ====================================================" << G4endl;
  G4cout << "   Beam Particle: " 
         << fParam->GetBeamParticle()->GetParticleName() << G4endl
         << "   Ekin(MeV)    = " << fParam->GetBeamEnergy()/MeV
         << G4endl
         << "   Z(mm)        = " << fParam->GetPositionZ()/mm 
         << G4endl;
  G4cout << " ================== run summary =====================" << G4endl;
  G4int prec = G4cout.precision(5);
  G4cout << "   End of Run TotNbofEvents    = " 
         << nEvt << G4endl;
  G4cout << "   Energy(keV) per ADC channel = " 
         << 1.0/(keV*fFactorALICE) << G4endl;

  G4cout << G4endl;
  G4cout << "   Mean energy deposit in absorber = " <<
    y1/keV << " +- " << y2*std::sqrt(norm)/keV << " keV; ";
  if(y1 > 0.0) { G4cout << "   RMS/Emean = " << y2/y1; }
  G4cout << G4endl;
  G4cout << "   Mean number of steps in absorber= " 
         << fTotStepGas << ";  mean number of ion-clusters = " 
         << fTotCluster << " MeanCluster= " << fMeanCluster    
         << G4endl;
  G4cout << G4endl;

  G4cout << " ====== Energy deposit distribution   Noverflows= " << fOverflow 
         << " ====== " << G4endl ;
  G4cout << " bin nb      Elow      entries     normalized " << G4endl;

  std::ofstream fileOut("distribution.out", std::ios::out );
  fileOut.setf( std::ios::scientific, std::ios::floatfield );

  x1 = 0.0;

  fileOut << fNbins << G4endl;
 
  for(G4int j=0; j<fNbins; ++j) 
  {
    G4cout << std::setw(5) << j << std::setw(10) << x1/keV 
           << std::setw(12) << fEgas[j] << std::setw(12) << fEgas[j]*norm 
           << G4endl ;
    fileOut << x1/keV << "\t" << fEgas[j] << G4endl;
    x1 += de;
  }
  G4cout.precision(prec);
 
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // normalize histograms
  G4double normf = fParam->GetNormFactor();
  analysisManager->ScaleH1(1,norm);
  analysisManager->ScaleH1(2,norm);
  analysisManager->ScaleH1(3,norm*normf);
 
  G4cout << " ================== run end ==========================" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::BeginOfEvent()
{
  fTotEdep = 0.0;
  fStepGas = 0;
  fCluster = 0;
  ++fEvt;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfEvent()
{
  fTotStepGas += fStepGas;
  fTotCluster += fCluster;

  if(fWidthALICE > 0.0) {
    G4double x = G4RandGauss::shoot(0.,fWidthALICE);
    fTotEdep += x;
    fTotEdep = std::max(fTotEdep, 0.0);
  }

  G4int idx = G4int(fTotEdep*fNbins/fMaxEnergy);

  if(idx < 0) { fEgas[0] += 1.0; }
  if(idx >= fNbins) { fOverflow += 1.0; }
  else { fEgas[idx] += 1.0; }
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // fill histo
  analysisManager->FillH1(1,fTotEdep/keV,1.0);
  analysisManager->FillH1(2,fCluster,1.0);
  analysisManager->FillH1(3,fTotEdep*fFactorALICE,1.0);
  fEdep.fill(fTotEdep, 1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);

  fTotStepGas  += localRun->fTotStepGas;
  fTotCluster  += localRun->fTotCluster;
  fMeanCluster += localRun->fMeanCluster;
  fOverflow    += localRun->fOverflow;

  G4StatDouble* stat = const_cast<G4StatDouble*>(localRun->GetStat());

  fEdep.add(stat);
 
  for(G4int j=0; j<fNbins; ++j)
  {
    fEgas[j] += localRun->fEgas[j]; 
  }
  
  G4Run::Merge(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddEnergy(G4double edep, const G4Step* step)
{
  if(1 < fVerbose) {
    G4cout << "Run::AddEnergy: e(keV)= " << edep/keV
           << G4endl;
  }
  fTotEdep += edep;
  if(step) {
    if(1 == step->GetTrack()->GetTrackID()) { fStepGas += 1.0; }

    fMeanCluster += fElIonPair->MeanNumberOfIonsAlongStep(step);
    fCluster += fElIonPair->SampleNumberOfIonsAlongStep(step);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
