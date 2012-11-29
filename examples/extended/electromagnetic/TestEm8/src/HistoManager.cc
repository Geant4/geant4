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
/// \file electromagnetic/TestEm8/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   HistoManager
//
// Author:      V.Ivanchenko 01.09.2010
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "Histo.hh"
#include "G4Step.hh"
#include "G4LossTableManager.hh"
#include "G4ElectronIonPair.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager* HistoManager::fManager = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager* HistoManager::GetPointer()
{
  if(!fManager) {
    fManager = new HistoManager();
  }
  return fManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager::HistoManager()
{
  fNHisto      = 2;
  fVerbose     = 1;
  fMaxEnergy   = 100.*keV;
  fBinsE       = 100;
  fBinsCluster = 200;

  fHisto     = new Histo();
  bookHisto();
  fElIonPair = G4LossTableManager::Instance()->ElectronIonPair();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager::~HistoManager()
{
  delete fHisto;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::bookHisto()
{
  fHisto->Add1D("10","Energy deposition in detector (keV)",
                fBinsE,0.0,fMaxEnergy/keV,1.0);
  fHisto->Add1D("11","Number of primary clusters",
                fBinsCluster,-0.5,fBinsCluster-0.5,1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BeginOfRun()
{
  // initilise scoring
  fEvt = 0.0;
  fTotStepGas = 0.0;
  fTotCluster = 0.0;
  fOverflow = 0.0;

  if(fHisto->IsActive()) { 
    fHisto->Activate(0, true); 
    fHisto->Activate(1, true); 
  }
  fHisto->Book();

  fEgas.resize(fBinsE,0.0);

  if(fVerbose > 0) {
    G4cout << "HistoManager: Histograms are booked and run has been started"
           << G4endl;
    G4cout << " BinsCluster= " << fBinsCluster << "    BinsE= " <<  fBinsE
           << "   Emax(keV)= " << fMaxEnergy/keV << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfRun()
{
  G4double norm = fEvt;
  if(fEvt > 0.0) { norm = 1.0/norm; }

  fTotStepGas *= norm;
  fTotCluster *= norm;
  fOverflow   *= norm;

  G4double y1 = 0.0;
  G4double y2 = 0.0;

  G4double de = fMaxEnergy/G4double(fBinsE);  
  G4double x1 = -de*0.5; 

  for(G4int j=0; j<fBinsE; ++j) {
    x1 += de;
    G4double y = fEgas[j]*norm;
    y1 += y*x1;
    y2 += y*x1*x1;
  }
  G4double norm1 = G4double(fBinsE);
  if(norm1 > 0.0) { norm1 = 1.0/norm1; }
  y1 *= norm1;
  y2 *= norm1;
  y2 -= y1*y1;
  if(y2 >= 0.0) { y2 = std::sqrt(y2); }

  G4cout << " ================== run summary =====================" << G4endl;
  G4int prec = G4cout.precision(5);
  G4cout << "   End of Run TotNbofEvents = " <<
           G4int(fEvt) << G4endl ;
  G4cout << G4endl;
  G4cout << "   Mean energy deposit in absorber = " <<
           y1/keV << " +- " << y2*std::sqrt(norm)/keV << " keV; ";
  if(y1 > 0.0) { G4cout << "   RMS/Emean = " << y2/y1; }
  G4cout << G4endl;
  G4cout << "   Mean number of steps in absorber= " 
         << fTotStepGas << ";  mean number of ion-clusters = " << fTotCluster    
         << G4endl;
  G4cout << G4endl;

  G4cout << " ====== Energy deposit distribution   Noverflows= " << fOverflow 
         << " ====== " << G4endl ;
  G4cout << " bin nb      Elow      entries     normalized " << G4endl;

  std::ofstream fileOut("distribution.out", std::ios::out );
  fileOut.setf( std::ios::scientific, std::ios::floatfield );

  x1 = 0.0;

  fileOut << fBinsE << G4endl;
 
  for(G4int j=0; j<fBinsE; ++j) 
  {
    G4cout << std::setw(5) << j << std::setw(10) << x1/keV 
           << std::setw(12) << fEgas[j] << std::setw(12) << fEgas[j]*norm 
           << G4endl ;
    fileOut << x1/keV << "\t" << fEgas[j] << G4endl;
    x1 += de;
  }
  G4cout.precision(prec);
  // normalise histograms
  fHisto->ScaleH1(0,norm);
  fHisto->ScaleH1(1,norm);
  
  fHisto->Save();
  G4cout << " ================== run end ==========================" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BeginOfEvent()
{
  fEvt += 1.0;
  fTotEdep = 0.0;
  fStepGas = 0;
  fCluster = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfEvent()
{
  fTotStepGas += fStepGas;
  fTotCluster += fCluster;
  
  G4int idx = G4int(fTotEdep*fBinsE/fMaxEnergy);
  if(idx < 0) { idx = 0; }
  if(idx >= fBinsE) { fOverflow += 1.0; }
  else { fEgas[idx] += 1.0; }

  // fill histo
  fHisto->Fill(0,fTotEdep/keV,1.0);
  fHisto->Fill(1,fCluster,1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::AddEnergy(G4double edep, G4Step* step)
{
  if(1 < fVerbose) {
    G4cout << "HistoManager::AddEnergy: e(keV)= " << edep/keV
           << G4endl;
  }
  fTotEdep += edep;
  if(step->GetTrack()->GetTrackID() == 1) { fStepGas += 1.0; }

  fCluster += fElIonPair->SampleNumberOfIonsAlongStep(step);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

