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
/// \file electromagnetic/TestEm5/src/Run.cc
/// \brief Implementation of the Run class
//
// $Id: Run.cc 71376 2013-06-14 07:44:50Z maire $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"

#include "G4EmCalculator.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det, PrimaryGeneratorAction* prim, 
         HistoManager* histoMgr)
:fDetector(det), fPrimary(prim), fHistoManager(histoMgr)
{
  fAnalysisManager = G4AnalysisManager::Instance();

  G4double length  = fDetector->GetAbsorSizeX();
  fOffsetX = -0.5 * length;

  fVerboseLevel = 1;
  fNevt = 0;
  fProjRange = fProjRange2 = 0.;

  fHistoManager->Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{ 

}

 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);

  fNevt   += localRun->GetNumberOfEvent();
  fProjRange += localRun->fProjRange;
  fProjRange2 += localRun->fProjRange2;
  
  G4Run::Merge(run); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun(double binLength)
{

#ifndef G4MULTITHREADED
  fNevt += this->GetNumberOfEvent();
#endif

  G4int nEvents = fNevt;
  if (nEvents == 0) { return; }

  //run conditions
  //  
  const G4Material* material = fDetector->GetAbsorMaterial();
  G4double density = material->GetDensity();
  G4String matName = material->GetName();
  const G4ParticleDefinition* part = 
    fPrimary->GetParticleGun()->GetParticleDefinition();
  G4String particle = part->GetParticleName();    
  G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();

  if(GetVerbose() > 0){
    G4cout << "\n The run consists of " << nEvents << " "<< particle << " of "
           << G4BestUnit(energy,"Energy") << " through " 
           << G4BestUnit(fDetector->GetAbsorSizeX(),"Length") << " of "
           << matName << " (density: " 
           << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
    G4cout<<"Proj "<<fProjRange<<" "<<fProjRange2<<G4endl;
  };
         
  //compute projected range and straggling

  fProjRange /= nEvents; fProjRange2 /= nEvents;
  G4double rms = fProjRange2 - fProjRange*fProjRange;        
  if (rms>0.) { rms = std::sqrt(rms); } 
  else { rms = 0.; }

  if(GetVerbose() > 0){
    G4cout.precision(5);       
    G4cout << "\n projected Range= " << G4BestUnit(fProjRange, "Length")
           << "   rms= "             << G4BestUnit(rms, "Length")
           << G4endl;
  };

  G4double ekin[100], dedxproton[100], dedxmp[100];
  G4EmCalculator calc;
  calc.SetVerbose(2);
  G4int i;
  for(i = 0; i < 100; ++i) {
    ekin[i] = std::pow(10., 0.1*G4double(i)) * keV;
    dedxproton[i] = 
      calc.ComputeElectronicDEDX(ekin[i], "proton", matName);
    dedxmp[i] = 
      calc.ComputeElectronicDEDX(ekin[i], "monopole", matName);
  }

  if(GetVerbose() > 0){
    G4cout << "### Stopping Powers" << G4endl;
    for(i=0; i<100; i++) {
      G4cout << " E(MeV)= " << ekin[i] << "  dedxp= " << dedxproton[i]
             << " dedxmp= " << dedxmp[i]
             << G4endl;
    }
  }
  G4cout << "### End of stopping power table" << G4endl;

  // normalize histogram
  G4double fac = (mm/MeV) / (nEvents * binLength);
  fHistoManager->Scale(1,fac);

  if(GetVerbose() > 0){
    G4cout << "Range table for " << matName << G4endl;
  }

  for(i=0; i<100; ++i) {
    G4double e = std::log10(ekin[i] / MeV) + 0.05;
    fHistoManager->FillHisto(2, e, dedxproton[i]);
    fHistoManager->FillHisto(3, e, dedxmp[i]);
    fHistoManager->FillHisto(4, e, 
                   std::log10(calc.GetRange(ekin[i],"proton",matName)/mm));
    fHistoManager->FillHisto(5, e, 
                   std::log10(calc.GetRange(ekin[i],"monopole",matName)/mm));
  }
  

  if(fAnalysisManager) {

    if(fAnalysisManager->IsActive()) {

      // Write histogram file
      if(!fAnalysisManager->Write()) {
        G4Exception ("Histo::Save()", "hist01", FatalException, 
                     "Cannot write ROOT file.");
      }
      G4cout << "### Histo::Save: Histograms are saved" << G4endl;
      if(fAnalysisManager->CloseFile() && fVerboseLevel) {
        G4cout << "                 File is closed" << G4endl;
      }
    }
    
    delete fAnalysisManager;
    fAnalysisManager = 0;
  }


}   

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void Run::FillHisto(G4int histoId, G4double v1, G4double v2)
{
  if(fAnalysisManager) 
    fHistoManager->FillHisto(histoId, v1, v2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
