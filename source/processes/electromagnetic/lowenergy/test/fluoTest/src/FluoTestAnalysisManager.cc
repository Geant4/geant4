#ifdef  G4ANALYSIS_USE

#include <stdlib.h>
#include "g4std/fstream"
#include "FluoTestAnalysisManager.hh"
#include "G4VAnalysisSystem.hh"
#include "FluoTestDetectorConstruction.hh"
#include "FluoTestAnalysisMessenger.hh"
#include "G4ios.hh"
#include "Interfaces/IHistoManager.h"
#include "Interfaces/IHistogram1D.h"
#include "Interfaces/IHistogram2D.h"

#include "NtupleTag/LizardNTupleFactory.h"
using namespace Lizard;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestAnalysisManager::FluoTestAnalysisManager(FluoTestDetectorConstruction* FluoTestDC): 
 
  histoManager(0),
  factory(0), 
  ntuple(0),
  Detector(FluoTestDC),
  histoGamDet(0),
  histoGamDetPre(0),
  histoGamLeavSam(0),
  histoEleLeavSam(0),
  //histoGamLS(0),
  //histoGamLSP(0),
  histoGamBornSam(0),
  histoEleBornSam(0),
  histoProtLeavSam(0),
  histoProtDetPre(0),
  histoSpectrum(0)
{
   // Define the messenger and the analysis system
  analysisMessenger = new FluoTestAnalysisMessenger(this);
 
  histoManager = createIHistoManager(); 
  assert (histoManager != 0);
  histoManager->selectStore("fluoTestHisto.hbook");
  factory = createNTupleFactory();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestAnalysisManager::~FluoTestAnalysisManager() 
{

  delete analysisMessenger; 
  delete  histoManager;
  delete ntuple;
  delete factory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestAnalysisManager::InsGamDet(double gD)
{
  histoGamDet->fill(gD);
}

void FluoTestAnalysisManager::InsGamBornSample(double gBs)
{
  histoGamBornSam->fill(gBs);
}
void FluoTestAnalysisManager::InsEleBornSample(double eBs)
{
  histoEleBornSam->fill(eBs);
}
/*
void FluoTestAnalysisManager::InsGamLS(double gS)
{
  histoGamLS->fill(gS);
}
void FluoTestAnalysisManager::InsGamLSP(double gSP)
{
  histoGamLSP->fill(gSP);
}*/
void FluoTestAnalysisManager::InsGamLeavSam(double gLs)
{
  histoGamLeavSam->fill(gLs);
}
void FluoTestAnalysisManager::InsEleLeavSam(double eLs)
{
  histoEleLeavSam->fill(eLs);
}

void FluoTestAnalysisManager::InsGamDetPre(double gDpr)
{
  histoGamDetPre->fill(gDpr);
}

void FluoTestAnalysisManager::InsSpectrum(double sp)
{
  histoSpectrum->fill(sp);
}
void FluoTestAnalysisManager::InsProtDetPre(double pdp)
{
  histoProtDetPre->fill(pdp);
}
void FluoTestAnalysisManager::InsProtLeavSam(double pLs)
{
  histoProtLeavSam->fill(pLs);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// This member reset the histograms and it is called at the begin
// of each run; here we put the inizialization so that the histograms have 
// always the right dimensions depending from the detector geometry

void FluoTestAnalysisManager::BeginOfRun() 
{ 
  ntuple = factory->createC( "fluoTestHisto1.hbook::1" );
  assert ( ntuple != 0 );
  histoGamDet = histoManager->create1D("10","Energy deposit in the detector", 1000,0,6.5);
  histoGamDetPre = histoManager->create1D("20","Gammas reaching the detector", 1000,0,6.5);
  histoGamLeavSam = histoManager->create1D("30","Gammas leaving the sample", 1000,0,6.5);
  histoEleLeavSam = histoManager->create1D("40","Electrons leaving the sample", 1000,0,6.5);
  // histoGamLS = histoManager->create1D("50","Theta of gammas leaving the sample", 800,0.,pi);
  // histoGamLSP = histoManager->create1D("60","Phi of gammas leaving the sample", 800,-pi,pi);
  histoGamBornSam = histoManager->create1D("70"," Gammas born in the sample", 1000,0,6.5);
  histoEleBornSam = histoManager->create1D("80"," Electrons born in the sample", 1000,0,6.5);
  histoSpectrum =  histoManager->create1D("90","Spectrum of the incident photons",1000,0,6.5);
  histoProtLeavSam = histoManager->create1D("50","Gammas leaving the sample", 1000,0,6.5);
  histoProtDetPre = histoManager->create1D("60","Gammas reaching the detector", 1000,0,6.5);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
//   This member is called at the end of each run 

void FluoTestAnalysisManager::EndOfRun(G4int n) 
{
  histoManager->store("10");
  histoManager->store("20");
  histoManager->store("30");
  histoManager->store("40");
  histoManager->store("50");
  histoManager->store("60");
  histoManager->store("70");
  histoManager->store("80");
  histoManager->store("90");


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//  This member is called at the end of every event 
void FluoTestAnalysisManager::EndOfEvent(G4int flag) 
{
  
}

#endif







