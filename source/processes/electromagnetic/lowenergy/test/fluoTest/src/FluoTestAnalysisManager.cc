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
  //histoDetETot(0),
  histoGamDetPre(0),
  //histoGamDetPost(0),
  histoGamLeavSam(0),
  histoEleLeavSam(0),
    histoGamLS(0),
  histoGamLSP(0),
  histoGamBornSam(0),
  histoEleBornSam(0)
  // histoOtherPartDet(0),
 
{
  // Define the messenger and the analysis system
  analysisMessenger = new FluoTestAnalysisMessenger(this);
 
  histoManager = createIHistoManager(); 
  assert (histoManager != 0);
  histoManager->selectStore("fluoTestHisto.hbook");

  factory = createNTupleFactory();
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestAnalysisManager::~FluoTestAnalysisManager() {
 
  delete analysisMessenger; 
   delete  histoManager;
   delete ntuple;
  delete factory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  void FluoTestAnalysisManager::InsGamDet(G4double gD)
  {
  histoGamDet->fill(gD);
  }
/* 
  void FluoTestAnalysisManager::InsDetETot(G4double dEt)
  {
  histoDetETot->fill(dEt);
  }
*/
void FluoTestAnalysisManager::InsGamBornSample(G4double gBs)
{
  histoGamBornSam->fill(gBs);
}
void FluoTestAnalysisManager::InsEleBornSample(G4double eBs)
{
  histoEleBornSam->fill(eBs);
}

void FluoTestAnalysisManager::InsGamLS(G4double gS)
{
  histoGamLS->fill(gS);
}
void FluoTestAnalysisManager::InsGamLSP(G4double gSP)
{
  histoGamLSP->fill(gSP);
}
void FluoTestAnalysisManager::InsGamLeavSam(G4double gLs)
{
  histoGamLeavSam->fill(gLs);
}
void FluoTestAnalysisManager::InsEleLeavSam(G4double eLs)
{
  histoEleLeavSam->fill(eLs);
}

  void FluoTestAnalysisManager::InsGamDetPre(G4double gDpr)
  {
   histoGamDetPre->fill(gDpr);
   }
/*

  void FluoTestAnalysisManager::InsGamDetPost(G4double gDps)
  {
  histoGamDetPost->fill(gDps);
  }
  
void FluoTestAnalysisManager::InsOtherPart(G4double oP)
{
histoOtherPartDet->fill(oP);
}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// This member reset the histograms and it is called at the begin
// of each run; here we put the inizialization so that the histograms have 
// always the right dimensions depending from the detector geometry

void FluoTestAnalysisManager::BeginOfRun() 
{
  ntuple = factory->createC( "fluoTestHisto1.hbook::1" );
  assert ( ntuple != 0 );
  histoGamDet = histoManager->create1D("10","Energy deposit in the detector", 100,0.,10.);
  histoGamDetPre = histoManager->create1D("20","Gammas reaching the detector", 100,0.,10.);
  histoGamLeavSam = histoManager->create1D("30","Gammas leaving the sample", 100,0.,10.);
  histoEleLeavSam = histoManager->create1D("40","Electrons leaving the sample", 100,0.,10.);
  histoGamLS = histoManager->create1D("50","Theta of gammas leaving the sample", 100,0.,pi);
  histoGamLSP = histoManager->create1D("60","Phi of gammas leaving the sample", 100,-pi,pi);
  histoGamBornSam = histoManager->create1D("70"," Gammas born in the sample", 100,0.,10.);
  histoEleBornSam = histoManager->create1D("80"," Electrons born in the sample", 100,0.,10.);


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

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//  This member is called at the end of every event 
void FluoTestAnalysisManager::EndOfEvent(G4int flag) 
{

}

#endif







