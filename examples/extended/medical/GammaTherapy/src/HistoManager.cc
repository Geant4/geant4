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
/// \file electromagnetic/TestEm5/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
// $Id: HistoManager.cc 102694 2017-02-17 09:07:55Z maire $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4SystemOfUnits.hh"
#include "DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
{

  fNBinsZ    = 60;
  fNBinsR    = 80;
  fNBinsE    = 200;

  fAbsorberZ = 300.*mm;
  fAbsorberR = 200.*mm;
  fScoreZ    = 100.*mm;
  fMaxEnergy = 50.*MeV;

  fStepZ = fStepR = fStepE = 0.0;

  Book();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);    // enable inactivation of histograms

  // Creating an 1-dimensional histograms in the root directory of the tree
  fHisto.assign(10,0);
  int iHisto=0;
  fHisto[iHisto] = analysisManager->CreateH1("10",
    "Energy deposit at radius (mm) normalised on 1st channel",
    fNBinsR, 0., fAbsorberR/mm);
  
  iHisto++;
  fHisto[iHisto] = analysisManager->CreateH1("11",
    "Energy deposit at radius (mm) normalised to integral",
    fNBinsR, 0., fAbsorberR/mm);

  iHisto++;
  fHisto[iHisto] = analysisManager->CreateH1("12",
    "Energy deposit (MeV/kg/electron) at radius (mm)",
    fNBinsR, 0., fAbsorberR/mm);

  iHisto++;
  fHisto[iHisto] = analysisManager->CreateH1("13",
    "Energy profile (MeV/kg/electron) over Z (mm)",fNBinsZ,0.,fAbsorberZ/mm);

  iHisto++;
  fHisto[iHisto] = analysisManager->CreateH1("14",
    "Energy profile (MeV/kg/electron) over Z (mm) at Central Voxel",
    fNBinsZ, 0., fAbsorberZ/mm);

  iHisto++;
  fHisto[iHisto] = analysisManager->CreateH1("15",
    "Energy (MeV) of fGamma produced in the target",
    fNBinsE, 0., fMaxEnergy/MeV);

  iHisto++;
  fHisto[iHisto] = analysisManager->CreateH1("16",
    "Energy (MeV) of fGamma before phantom",fNBinsE,0.,fMaxEnergy/MeV);

  iHisto++;
  fHisto[iHisto] = analysisManager->CreateH1("17",
    "Energy (MeV) of electrons produced in phantom",fNBinsE,0.,fMaxEnergy/MeV);

  iHisto++;
  fHisto[iHisto] = analysisManager->CreateH1("18",
    "Energy (MeV) of electrons produced in target",fNBinsE,0.,fMaxEnergy/MeV);

  iHisto++;
  fHisto[iHisto] = analysisManager->CreateH1("19",
    "Gamma Energy Fluence (MeV/cm2) at radius(mm) in front of phantom",
    fNBinsR, 0., fAbsorberR/mm);

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for(int i=0; i<iHisto+1; i++) analysisManager->SetH1Activation(i, false);

}


void HistoManager::Update(DetectorConstruction* det, bool bForceActivation)
{

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  if(bForceActivation) {
    for(int i=0; i<(int)fHisto.size(); i++) 
      analysisManager->SetH1Activation(fHisto[i], true);
    analysisManager->SetActivation(true);
  }

  if(analysisManager->IsActive()) {

    // Check nBinsR / fAbsorberR histograms
    if(det->GetNumberDivR()!=fNBinsR||
       fabs(det->GetAbsorberR()-fAbsorberR)>0.01*mm)
      {
        fNBinsR = det->GetNumberDivR();
        fAbsorberR = det->GetAbsorberR();
        std::vector<G4int> histoId { 0, 1, 2, 9 };
        for(auto v : histoId)
          {
            analysisManager->SetH1(fHisto[v], fNBinsR, 0., fAbsorberR/mm);
          }

      }
    
    // Check nBinsZ / fAbsorberZ histograms
    if(det->GetNumberDivZ()!=fNBinsZ||
       fabs(det->GetAbsorberZ()-fAbsorberZ)>0.01*mm)
      {
        fNBinsZ = det->GetNumberDivZ();
        fAbsorberZ = det->GetAbsorberZ();
        std::vector<G4int> histoId { 3, 4 };
        for(auto v : histoId)
          {
            analysisManager->SetH1(fHisto[v], fNBinsZ, 0., fAbsorberZ/mm);
          }
      }
    
    // Check nBinsE / fAbsorberE histograms
    if(det->GetNumberDivE()!=fNBinsE||
       fabs(det->GetMaxEnergy()-fMaxEnergy)>0.01)
      {
        fNBinsE = det->GetNumberDivE();
        fMaxEnergy = det->GetMaxEnergy();
        std::vector<G4int> histoId { 5, 6, 7 ,8 };
        for(auto v : histoId)
          {
            analysisManager->SetH1(fHisto[v], fNBinsE, 0., fMaxEnergy/MeV);
          }
      }
  }

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::DumpHistoParameters()
{
  if(!G4AnalysisManager::Instance()->IsActive()) return;

  for(int i=0; i<(int)fHisto.size(); i++)
    {
      G4int histoId=fHisto[i];      
      G4String title = G4AnalysisManager::Instance()->GetH1Title(histoId);
      G4int nbins = G4AnalysisManager::Instance()->GetH1Nbins(histoId);
      G4double xmin = G4AnalysisManager::Instance()->GetH1Xmin(histoId);
      G4double xmax = G4AnalysisManager::Instance()->GetH1Xmax(histoId);
      G4double width = G4AnalysisManager::Instance()->GetH1Width(histoId);
      G4cout<<"Histogram parameters : "<<i<<" "<<histoId<<" : "<<nbins<<" ";
      G4cout<<xmin<<"/"<<xmax<<" "<<width<<G4endl;
      
    }

}




