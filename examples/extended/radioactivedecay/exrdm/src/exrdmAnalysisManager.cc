//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "exrdmAnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "exrdmHisto.hh"

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
#ifdef G4ANALYSIS_USE
   bookHisto();
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exrdmAnalysisManager::~exrdmAnalysisManager()
{
  delete histo;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmAnalysisManager::bookHisto()
{
  histEMax = 15.0*MeV;
  histEMin = .0*MeV;
  histNBin = 100;

  histo->add1D("10",
    "Energy deposit (MeV) in the traget",histNBin,histEMin,histEMax,MeV);
  histo->add1D("11",
    "Energy deposit (MeV) in the detector",histNBin,histEMin,histEMax,MeV);
  histo->add1D("12",
    "Total energy spectrum (MeV) of the traget and detector",histNBin,histEMin,histEMax,MeV);
  histo->add1D("13",
    "Coincidence spectrum (MeV) between the traget and detector",histNBin,histEMin,histEMax,MeV);
  histo->add1D("14",
    "Anti-coincidence spectrum (MeV) in the traget",histNBin,histEMin,histEMax,MeV);
  histo->add1D("15",
    "Anti-coincidence spectrum (MeV) in the detector",histNBin,histEMin,histEMax,MeV);
  histo->add1D("16",
	       "Decay emission spectrum (MeV)",histNBin,histEMin,histEMax,MeV);
  // in aida these histos are indiced from 0-6
  //
  histo->addTuple( "100", "Emitted Particles","string Name, float Energy, Time, Weight" );
  histo->addTuple( "200", "RadioIsotopes","string Name, float Time, Weight" );
  histo->addTuple( "300", "Energy Depositions","float Energy, Time, Weight" );

  // at the moment anaphe can only handle float ntuple
  //histo->addTuple( "100", "Emitted Particles"," float Name, Energy, Time, Weight" );
  //histo->addTuple( "200", "RadioIsotopes"," float Name, Energy, Time, Weight" );
  //histo->addTuple( "300", "Energy Depositions"," float Energy, Time, Weight" );
  //    the ntuples are indeced from 0-2 in aida
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmAnalysisManager::BeginOfRun()
{
#ifdef G4ANALYSIS_USE
  histo->book();
#endif
  if(verbose > 0) {
    G4cout << "exrdmAnalysisManager: Histograms are booked and the run has been started"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmAnalysisManager::EndOfRun()
{
#ifdef G4ANALYSIS_USE
  histo->save();  
#endif
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
  histo->fillTuple(2,"Energy", edep/MeV);
  histo->fillTuple(2,"Weight",weight);
  histo->fillTuple(2,"Time",time/second);
  histo->addRow(2);
  // 
  exrdmEnergyDeposition A(edep,time,weight);
  Edepo.push_back(A);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmAnalysisManager::AddParticle(G4String particleName, G4double energy, G4double weight, G4double time )
{
  if(1 < verbose) {
    G4cout << "exrdmAnalysisManager::AddParticle: " << particleName
           << G4endl;
  }
  histo->fillTuple(0,"Name", particleName);
  histo->fillTuple(0,"Energy",energy/MeV);
  histo->fillTuple(0,"Weight",weight);
  histo->fillTuple(0,"Time",time/second);
  histo->addRow(0);
  // now fill th emission spectrum
  if (energy) histo->fillHisto(6,energy/MeV,weight);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmAnalysisManager::AddIsotope(G4String particleName,G4double weight, G4double time )
{
  if(1 < verbose) {
    G4cout << "exrdmAnalysisManager::AddIsotope: " << particleName
           << G4endl;
  }
  histo->fillTuple(1,"Name",particleName);
  histo->fillTuple(1,"Weight",weight);
  histo->fillTuple(1,"Time",time/second);
  histo->addRow(1);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

