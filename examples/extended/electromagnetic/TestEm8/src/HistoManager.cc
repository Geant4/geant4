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
// $Id: HistoManager.cc,v 1.4 2010-11-19 10:15:33 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  nHisto       = 2;
  verbose      = 1;
  maxEnergy    = 100.*keV;
  nBinsE       = 100;
  nBinsCluster = 200;

  histo   = new Histo();
  elIonPair = G4LossTableManager::Instance()->ElectronIonPair();

  bookHisto();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager::~HistoManager()
{
  delete histo;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::bookHisto()
{
  histo->add1D("10",
    "Energy deposition in detector (keV)",nBinsE,0.0,maxEnergy/keV,1.0);

  histo->add1D("11",
    "Number of primary clusters",nBinsCluster,-0.5,G4double(nBinsCluster)-0.5,1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BeginOfRun()
{
  // initilise scoring
  nEvt = 0;
  nTotStepGas = 0;
  nTotCluster = 0;
  overflow = 0;

  histo->book();

  Egas.resize(nBinsE,0.0);

  if(verbose > 0) {
    G4cout << "HistoManager: Histograms are booked and run has been started"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfRun()
{
  G4double norm = nEvt;
  if(nEvt > 0.0) { norm = 1.0/norm; }

  nTotStepGas *= norm;
  nTotCluster *= norm;
  overflow    *= norm;

  G4double y1 = 0.0;
  G4double y2 = 0.0;

  G4double de = maxEnergy/G4double(nBinsE);  
  G4double x1 = -de*0.5; 

  for(G4int j=0; j<nBinsE; ++j) {
    x1 += de;
    G4double y = Egas[j]*norm;
    y1 += y*x1;
    y2 += y*x1*x1;
  }
  y2 -= y1*y1;
  if(y2 >= 0.0) { y2 = std::sqrt(y2); }

  G4cout << " ================== run summary =====================" << G4endl;
  G4int prec = G4cout.precision(5);
  G4cout << "   End of Run TotNbofEvents = " <<
           G4int(nEvt) << G4endl ;
  G4cout << G4endl;
  G4cout << "   Mean energy deposit in absorber = " <<
           y1/keV << " +- " << y2*std::sqrt(norm)/keV << " keV; ";
  if(y1 > 0.0) { G4cout << "   RMS/Emean = " << y2/y1; }
  G4cout << G4endl;
  G4cout << "   Mean number of steps in absorber= " 
	 << nTotStepGas << ";  mean number of ion-clusters = " << nTotCluster    
         << G4endl;
  G4cout << G4endl;

  G4cout << " ====== Energy deposit distribution   Noverflows= " << overflow 
	 << " ====== " << G4endl ;
  G4cout << " bin nb      Elow      entries     normalized " << G4endl;

  std::ofstream fileOut("distribution.out", std::ios::out );
  fileOut.setf( std::ios::scientific, std::ios::floatfield );

  x1 = 0.0;

  fileOut<<nBinsE<<G4endl;
 
  for(G4int j=0; j<nBinsE; ++j) 
  {
    G4cout << std::setw(5) << j << std::setw(10) << x1/keV 
	   << std::setw(12) << Egas[j] << std::setw(12) << Egas[j]*norm 
	   << G4endl ;
    fileOut << x1/keV << "\t" << Egas[j] << G4endl;
    x1 += de;
  }
  G4cout.precision(prec);
  // normalise histograms
  histo->scale(0,norm);
  histo->scale(1,norm);
  
  histo->save();
  G4cout << " ================== run end ==========================" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BeginOfEvent()
{
  nEvt += 1;
  totEdep = 0.0;
  nStepGas = 0;
  nCluster = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfEvent()
{
  nTotStepGas += nStepGas;
  nTotCluster += nCluster;
  
  G4int idx = G4int(totEdep*nBinsE/maxEnergy);
  if(idx < 0) { idx = 0; }
  if(idx >= nBinsE) { overflow += 1.0; }
  else { Egas[idx] += 1.0; }

  // fill histo
  histo->fill(0,totEdep/keV,1.0);
  histo->fill(1,nCluster,1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::AddEnergy(G4double edep, G4Step* step)
{
  if(1 < verbose) {
    G4cout << "HistoManager::AddEnergy: e(keV)= " << edep/keV
           << G4endl;
  }
  totEdep += edep;
  if(step->GetTrack()->GetTrackID() == 1) { nStepGas += 1.0; }

  nCluster += elIonPair->SampleNumberOfIonsAlongStep(step);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

