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
/// \file radioactivedecay/rdecay02/src/exrdmAnalysisManager.cc
/// \brief Implementation of the exrdmAnalysisManager class
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "exrdmAnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "exrdmHisto.hh"
#include "G4ProcessTable.hh"
#include "G4RadioactiveDecay.hh"
#include "G4TwoVector.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exrdmAnalysisManager* exrdmAnalysisManager::fManager = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exrdmAnalysisManager* exrdmAnalysisManager::GetInstance()
{
  if(!fManager) {
    fManager = new exrdmAnalysisManager();
  }
  return fManager;
}
void exrdmAnalysisManager::Dispose()
{
  delete fManager;
  fManager = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exrdmAnalysisManager::exrdmAnalysisManager()
: fVerbose(0), fNEvt1(-1), fNEvt2(-2),
  fHistEMax (15.0*MeV), fHistEMin (0.),fHistNBin(100),
  fTargetThresE(10*keV), fDetectorThresE(10*keV),fPulseWidth(1.*microsecond)
{
  fEdepo.clear();
  fHisto   = new exrdmHisto();
  BookHisto();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exrdmAnalysisManager::~exrdmAnalysisManager()
{
#ifdef G4ANALYSIS_USE 
#define HISTDELETE
#endif
#ifdef G4ANALYSIS_USE_ROOT
#define HISTDELETE
#endif
#ifdef HISTDELETE
  delete fHisto;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmAnalysisManager::BookHisto()
{
  fHistEMax = 15.0*MeV;
  fHistEMin = .0*MeV;
  fHistNBin = 100;

  fHisto->Add1D("H10",
    "Energy deposit (MeV) in the traget",fHistNBin,fHistEMin,fHistEMax,MeV);
  fHisto->Add1D("H11",
    "Energy deposit (MeV) in the detector",fHistNBin,fHistEMin,fHistEMax,MeV);
  fHisto->Add1D("H12",
    "Total energy spectrum (MeV) of the traget and detector",fHistNBin,
    fHistEMin,fHistEMax,MeV);
  fHisto->Add1D("H13",
    "Coincidence spectrum (MeV) between the traget and detector",fHistNBin,
    fHistEMin,fHistEMax,MeV);
  fHisto->Add1D("H14",
    "Anti-coincidence spectrum (MeV) in the traget",fHistNBin,
    fHistEMin,fHistEMax,MeV);
  fHisto->Add1D("H15",
    "Anti-coincidence spectrum (MeV) in the detector",fHistNBin,
    fHistEMin,fHistEMax,MeV);
  fHisto->Add1D("H16",
               "Decay emission spectrum (MeV)",fHistNBin,fHistEMin,fHistEMax,MeV);
  // in aida these histos are indiced from 0-6
  //
  fHisto->AddTuple( "T1", "Emitted Particles",
                                           "double PID, Energy, Time, Weight" );
  fHisto->AddTuple( "T2", "RadioIsotopes","double PID, Time, Weight" );
  fHisto->AddTuple( "T3", "Energy Depositions","double Energy, Time, Weight" );
  fHisto->AddTuple( "RDecayProducts", "All Products of RDecay",
                    "double PID, Z, A, Energy, Time, Weight" );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmAnalysisManager::BeginOfRun()
{
  fHisto->Book();
  G4cout << "exrdmAnalysisManager: Histograms are booked and the run has been started" << G4endl;
  G4ProcessTable *pTable = G4ProcessTable::GetProcessTable();
  G4RadioactiveDecay * rDecay = (G4RadioactiveDecay *)
     pTable->FindProcess("RadioactiveDecay", "GenericIon");
  if (rDecay != NULL) {
    if (!rDecay->IsAnalogueMonteCarlo()) {
      std::vector<G4RadioactivityTable*> theTables =
         rDecay->GetTheRadioactivityTables();
      for (size_t i = 0 ; i < theTables.size(); i++) {
          theTables[i]->GetTheMap()->clear();
          G4cout << " Clear the radioactivity map: 0, new size = "
                 << theTables[i]->GetTheMap()->size() << G4endl;
      }
    }  
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmAnalysisManager::EndOfRun(G4int nevent)
{
  fHisto->Save();
  G4cout << "exrdmAnalysisManager: Histograms have been saved!" << G4endl;

  // Output the induced radioactivities
  //   in VR mode only
  //
  G4ProcessTable *pTable = G4ProcessTable::GetProcessTable();
  G4RadioactiveDecay * rDecay = (G4RadioactiveDecay *)
    pTable->FindProcess("RadioactiveDecay", "GenericIon");
  if (rDecay != NULL) {
    if (!rDecay->IsAnalogueMonteCarlo()) {
      G4String fileName = fHisto->GetFileName() + ".activity" ;
      std::ofstream outfile (fileName, std::ios::out );
      std::vector<G4RadioactivityTable*> theTables =
         rDecay->GetTheRadioactivityTables();
      for (size_t i = 0 ; i < theTables.size(); i++) {
            G4double rate, error;
            outfile << "Radioactivities in decay window no. " << i << G4endl;
            outfile <<
                   "Z \tA \tE \tActivity (decays/window) \tError (decays/window) "
                    << G4endl;
            map<G4ThreeVector,G4TwoVector> *aMap = theTables[i]->GetTheMap();
            map<G4ThreeVector,G4TwoVector>::iterator iter;
            for(iter=aMap->begin(); iter != aMap->end(); iter++) {
              rate = iter->second.x()/nevent;
              error = std::sqrt(iter->second.y())/nevent;
              if ( rate < 0.) rate = 0.; // statically it can be < 0. but it's unphysical
              outfile << iter->first.x() <<"\t"<< iter->first.y() <<"\t"
                          << iter->first.z() << "\t" << rate <<"\t" << error << G4endl;
            }
            outfile << G4endl;
      }
      outfile.close();
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmAnalysisManager::BeginOfEvent()
{
  fEdepo.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmAnalysisManager::EndOfEvent()
{
  if (fEdepo.size()) {
    std::sort(fEdepo.begin(),fEdepo.end());
    G4double TarE = fEdepo[0].GetEnergy();
    G4double Time = fEdepo[0].GetTime();
    G4double TarW = fEdepo[0].GetEnergy()*fEdepo[0].GetWeight();
    G4double DetE = 0.;
    G4double DetW = 0.;
    G4double ComW = 0.;
    if (TarE< 0.) {
      DetE = -TarE;
      DetW = -TarW;
      TarE = 0.;
      TarW = 0.;
    }
    for (size_t i = 1; i < fEdepo.size(); i++) {
      if (std::fabs((fEdepo[i].GetTime()- Time)/second) <= fPulseWidth) {
        if ( fEdepo[i].GetEnergy() > 0. ) {
          TarE += fEdepo[i].GetEnergy();
          TarW += fEdepo[i].GetEnergy()*fEdepo[i].GetWeight();
        } else {
          DetE -= fEdepo[i].GetEnergy();
          DetW -= fEdepo[i].GetEnergy()*fEdepo[i].GetWeight();
        }
      } else {
        // tallying
        if (TarE || DetE) ComW = (TarW+DetW)/(TarE+DetE);
        if (TarE) TarW /= TarE;
        if (DetE) DetW /= DetE;
        //
        if (TarE) fHisto->FillHisto(0,TarE,TarW); // target histogram
        if (DetE) fHisto->FillHisto(1,DetE,DetW); // Detector histogram
        if (TarE+DetE)  fHisto->FillHisto(2,DetE+TarE,ComW); // Total histogram
        if (DetE >= fDetectorThresE && TarE >= fTargetThresE )
          fHisto->FillHisto(3,DetE,DetW); // coincidence histogram
        if (TarE >= fTargetThresE && DetE < fDetectorThresE) 
          fHisto->FillHisto(4,TarE,TarW); // target anti-coincidence histogram
        if (TarE < fTargetThresE && DetE >= fDetectorThresE) 
          fHisto->FillHisto(5,DetE,DetW); // detector anti-coincidence histogram
        //
        //reset
        TarE = fEdepo[i].GetEnergy();
        Time = fEdepo[i].GetTime();
        TarW = fEdepo[i].GetEnergy()*fEdepo[i].GetWeight();
        DetE = 0.;
        DetW = 0.;
        if (TarE < 0) {
          DetE = -TarE;
          DetW = -TarW;
          TarE = 0.;
          TarW = 0.;
        }
      }
    }
    //tally the last hit
    if (TarE || DetE) ComW = (TarW+DetW)/(TarE+DetE);
    if (TarE) TarW /= TarE;
    if (DetE) DetW /= DetE;
    //
    if (TarE) fHisto->FillHisto(0,TarE,TarW); // target histogram
    if (DetE) fHisto->FillHisto(1,DetE,DetW); // Detector histogram
    if (TarE+DetE)  fHisto->FillHisto(2,DetE+TarE,ComW); // Total histogram
    if (DetE >= fDetectorThresE && TarE >= fTargetThresE )
      fHisto->FillHisto(3,DetE,DetW); // coincidence histogram
    if (TarE >= fTargetThresE && DetE < fDetectorThresE) 
      fHisto->FillHisto(4,TarE,TarW); // target anti-coincidence histogram
    if (TarE < fTargetThresE && DetE >= fDetectorThresE) 
      fHisto->FillHisto(5,DetE,DetW); // detector anti-coincidence histogram
    // now add zero energy to separat events
    AddEnergy(0.,0.,0.);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exrdmAnalysisManager::AddEnergy(G4double edep, G4double weight,
                                                                                                                           G4double time)
{
  if(1 < fVerbose) {
    G4cout << "exrdmAnalysisManager::AddEnergy: e(keV)= " << edep/keV 
           << " weight = " << weight << " time (s) = " <<  time/second
           << G4endl;
  }
  fHisto->FillTuple(2, 0, edep/MeV);
  fHisto->FillTuple(2,1,time/second);
  fHisto->FillTuple(2,2,weight);
  fHisto->AddRow(2);
  // 
  exrdmEnergyDeposition A(edep,time,weight);
  fEdepo.push_back(A);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmAnalysisManager::AddParticle(G4double pid, G4double energy,
                                       G4double weight, G4double time )
{
  if(1 < fVerbose) {
    G4cout << "exrdmAnalysisManager::AddParticle: " << pid
           << G4endl;
  }
  fHisto->FillTuple(0,0, pid);
  fHisto->FillTuple(0,1,energy/MeV);
  fHisto->FillTuple(0,2,time/second);
  fHisto->FillTuple(0,3,weight);
  fHisto->AddRow(0);
  // now fill th emission spectrum
  if (energy>0.0) fHisto->FillHisto(6,energy/MeV,weight);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmAnalysisManager::AddIsotope(G4double pid,G4double weight,
                                      G4double time )
{
  if(1 < fVerbose) {
    G4cout << "exrdmAnalysisManager::AddIsotope: " << pid
           << G4endl;
  }
  fHisto->FillTuple(1,0,pid);
  fHisto->FillTuple(1,1,time/second);
  fHisto->FillTuple(1,2,weight);
  fHisto->AddRow(1);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmAnalysisManager::AddDecayProduct(G4double pid,G4int Z, G4int A,
                                           G4double energy, G4double time,
                                           G4double weight)
{
  fHisto->FillTuple(3,0,pid);
  fHisto->FillTuple(3,1,double(Z));
  fHisto->FillTuple(3,2,double(A));
  fHisto->FillTuple(3,3,energy);
  fHisto->FillTuple(3,4,time/s);
  fHisto->FillTuple(3,5,weight);
  fHisto->AddRow(3);
}
