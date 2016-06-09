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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "exrdmAnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "exrdmHisto.hh"
#include "G4ProcessTable.hh"
#include "G4RadioactiveDecay.hh"

#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exrdmAnalysisManager* exrdmAnalysisManager::fManager = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exrdmAnalysisManager* exrdmAnalysisManager::getInstance()
{
  if(!fManager) {
    fManager = new exrdmAnalysisManager();
  }
  return fManager;
}
void exrdmAnalysisManager::dispose()
{
  delete fManager;
  fManager = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exrdmAnalysisManager::exrdmAnalysisManager()
{
  verbose = 0;
  nEvt1   = -1;
  nEvt2   = -1;
  targetThresE = 10*keV;
  detectorThresE = 10*keV;
  pulseWidth = 1.*microsecond;
  histo   = new exrdmHisto();
  bookHisto();
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
  delete histo;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmAnalysisManager::bookHisto()
{
  histEMax = 15.0*MeV;
  histEMin = .0*MeV;
  histNBin = 100;

  histo->add1D("H10",
    "Energy deposit (MeV) in the traget",histNBin,histEMin,histEMax,MeV);
  histo->add1D("H11",
    "Energy deposit (MeV) in the detector",histNBin,histEMin,histEMax,MeV);
  histo->add1D("H12",
    "Total energy spectrum (MeV) of the traget and detector",histNBin,histEMin,histEMax,MeV);
  histo->add1D("H13",
    "Coincidence spectrum (MeV) between the traget and detector",histNBin,histEMin,histEMax,MeV);
  histo->add1D("H14",
    "Anti-coincidence spectrum (MeV) in the traget",histNBin,histEMin,histEMax,MeV);
  histo->add1D("H15",
    "Anti-coincidence spectrum (MeV) in the detector",histNBin,histEMin,histEMax,MeV);
  histo->add1D("H16",
	       "Decay emission spectrum (MeV)",histNBin,histEMin,histEMax,MeV);
  // in aida these histos are indiced from 0-6
  //
  histo->addTuple( "T1", "Emitted Particles","double PID, Energy, Time, Weight" );
  histo->addTuple( "T2", "RadioIsotopes","double PID, Time, Weight" );
  histo->addTuple( "T3", "Energy Depositions","double Energy, Time, Weight" );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmAnalysisManager::BeginOfRun()
{
  histo->book();
  G4cout << "exrdmAnalysisManager: Histograms are booked and the run has been started" << G4endl;
  G4ProcessTable *pTable = G4ProcessTable::GetProcessTable();
  G4RadioactiveDecay * rDecay = (G4RadioactiveDecay *) pTable->FindProcess("RadioactiveDecay", "GenericIon");
  if (rDecay != NULL) {
    if (!rDecay->IsAnalogueMonteCarlo()) {
      std::vector<G4RadioactivityTable*> theTables =  rDecay->GetTheRadioactivityTables();
      for (size_t i = 0 ; i < theTables.size(); i++) {
	theTables[i]->GetTheMap()->clear();
	G4cout << " Clear the radioactivity map: 0, new size = " << theTables[i]->GetTheMap()->size() << G4endl;  
      }
    }  
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmAnalysisManager::EndOfRun(G4int nevent)
{
  histo->save();
  G4cout << "exrdmAnalysisManager: Histograms have been saved!" << G4endl;

  // Output the induced radioactivities
  //   in VR mode only
  //
  G4ProcessTable *pTable = G4ProcessTable::GetProcessTable();
  G4RadioactiveDecay * rDecay = (G4RadioactiveDecay *) pTable->FindProcess("RadioactiveDecay", "GenericIon");
  if (rDecay != NULL) {
    if (!rDecay->IsAnalogueMonteCarlo()) {
      G4String fileName = histo->getFileName() + ".activity" ;
      std::ofstream outfile (fileName, std::ios::out );
      std::vector<G4RadioactivityTable*> theTables =  rDecay->GetTheRadioactivityTables();
      for (size_t i = 0 ; i < theTables.size(); i++) {
	G4double rate;
	outfile << "Radioactivities in decay window no. " << i << G4endl;
	outfile << "Z \tA \tE \tActivity (decays/window) "<< G4endl;
	map<G4ThreeVector,G4double> *aMap = theTables[i]->GetTheMap();
	map<G4ThreeVector,G4double>::iterator iter;
	for(iter=aMap->begin(); iter != aMap->end(); iter++) {
	  rate = iter->second;
	  if ( rate < 0.) rate = 0.; // statically it can be < 0. but it's unphysical
	  outfile << iter->first.x() <<"\t"<< iter->first.y() <<"\t"<< iter->first.z() << "\t" << rate/nevent << G4endl;
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
  Edepo.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmAnalysisManager::EndOfEvent()
{
  if (Edepo.size()) {
    std::sort(Edepo.begin(),Edepo.end());
    G4double TarE = Edepo[0].GetEnergy();
    G4double Time = Edepo[0].GetTime();
    G4double TarW = Edepo[0].GetEnergy()*Edepo[0].GetWeight();
    G4double DetE = 0.;
    G4double DetW = 0.;
    G4double ComW = 0.;
    if (TarE< 0.) {
      DetE = -TarE;
      DetW = -TarW;
      TarE = 0.;
      TarW = 0.;
    }
    for (size_t i = 1; i < Edepo.size(); i++) {
      if (std::fabs((Edepo[i].GetTime()- Time)/second) <= pulseWidth) {
	if ( Edepo[i].GetEnergy() > 0. ) {
	  TarE += Edepo[i].GetEnergy();
	  TarW += Edepo[i].GetEnergy()*Edepo[i].GetWeight();
	} else {
	  DetE -= Edepo[i].GetEnergy();
	  DetW -= Edepo[i].GetEnergy()*Edepo[i].GetWeight();
	}
      } else {
	// tallying
	if (TarE || DetE) ComW = (TarW+DetW)/(TarE+DetE);
	if (TarE) TarW /= TarE;
	if (DetE) DetW /= DetE;
	//
	if (TarE) histo->fillHisto(0,TarE,TarW); // target histogram
	if (DetE) histo->fillHisto(1,DetE,DetW); // Detector histogram
	if (TarE+DetE)  histo->fillHisto(2,DetE+TarE,ComW); // Total histogram
	if (DetE >= detectorThresE && TarE >= targetThresE )
	  histo->fillHisto(3,DetE,DetW); // coincidence histogram
	if (TarE >= targetThresE && DetE < detectorThresE) 
	  histo->fillHisto(4,TarE,TarW); // target anti-coincidence histogram
	if (TarE < targetThresE && DetE >= detectorThresE) 
	  histo->fillHisto(5,DetE,DetW); // detector anti-coincidence histogram
	//
	//reset
	TarE = Edepo[i].GetEnergy();
	Time = Edepo[i].GetTime();
	TarW = Edepo[i].GetEnergy()*Edepo[i].GetWeight();
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
    if (TarE) histo->fillHisto(0,TarE,TarW); // target histogram
    if (DetE) histo->fillHisto(1,DetE,DetW); // Detector histogram
    if (TarE+DetE)  histo->fillHisto(2,DetE+TarE,ComW); // Total histogram
    if (DetE >= detectorThresE && TarE >= targetThresE )
      histo->fillHisto(3,DetE,DetW); // coincidence histogram
    if (TarE >= targetThresE && DetE < detectorThresE) 
      histo->fillHisto(4,TarE,TarW); // target anti-coincidence histogram
    if (TarE < targetThresE && DetE >= detectorThresE) 
      histo->fillHisto(5,DetE,DetW); // detector anti-coincidence histogram
    // now add zero energy to separat events
    AddEnergy(0.,0.,0.);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exrdmAnalysisManager::AddEnergy(G4double edep, G4double weight, G4double time)
{
  if(1 < verbose) {
    G4cout << "exrdmAnalysisManager::AddEnergy: e(keV)= " << edep/keV 
	   << " weight = " << weight << " time (s) = " <<  time/second
           << G4endl;
  }
  histo->fillTuple(2, 0, edep/MeV);
  histo->fillTuple(2,1,time/second);
  histo->fillTuple(2,2,weight);
  histo->addRow(2);
  // 
  exrdmEnergyDeposition A(edep,time,weight);
  Edepo.push_back(A);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmAnalysisManager::AddParticle(G4double pid, G4double energy, G4double weight, G4double time )
{
  if(1 < verbose) {
    G4cout << "exrdmAnalysisManager::AddParticle: " << pid
           << G4endl;
  }
  histo->fillTuple(0,0, pid);
  histo->fillTuple(0,1,energy/MeV);
  histo->fillTuple(0,2,time/second);
  histo->fillTuple(0,3,weight);
  histo->addRow(0);
  // now fill th emission spectrum
  if (energy>0.0) histo->fillHisto(6,energy/MeV,weight);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmAnalysisManager::AddIsotope(G4double pid,G4double weight, G4double time )
{
  if(1 < verbose) {
    G4cout << "exrdmAnalysisManager::AddIsotope: " << pid
           << G4endl;
  }
  histo->fillTuple(1,0,pid);
  histo->fillTuple(1,1,time/second);
  histo->fillTuple(1,2,weight);
  histo->addRow(1);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

