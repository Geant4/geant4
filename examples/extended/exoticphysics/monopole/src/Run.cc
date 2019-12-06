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
/// \file exoticphysics/monopole/src/Run.cc
/// \brief Implementation of the Run class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4EmCalculator.hh"
#include "G4Proton.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det, PrimaryGeneratorAction* prim)
 :fDetector(det), fPrimary(prim)
{
  fAnalysisManager = G4AnalysisManager::Instance();

  G4double length  = fDetector->GetAbsorSizeX();
  fOffsetX = -0.5 * length;

  fVerboseLevel = 1;
  fNevt = 0;
  fProjRange = fProjRange2 = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{}

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

  if ( ! G4Threading::IsMultithreadedApplication() ) {
    fNevt += this->GetNumberOfEvent();
  }

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
  const G4ParticleDefinition* proton = G4Proton::Proton(); 

  G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();

  if(GetVerbose() > 0){
    G4cout << "\n The run consists of " << nEvents << " "<< particle << " of "
           << G4BestUnit(energy,"Energy") << "\n through " 
           << G4BestUnit(fDetector->GetAbsorSizeX(),"Length") << " of "
           << matName << " (density: " 
           << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
    //G4cout<<"Proj "<<fProjRange<<" "<<fProjRange2<<G4endl;
  };
         
  //compute projected range and straggling
  fProjRange /= nEvents; fProjRange2 /= nEvents;
  G4double rms = fProjRange2 - fProjRange*fProjRange;        
  if (rms>0.) { rms = std::sqrt(rms); } 
  else { rms = 0.; }

  if(GetVerbose() > 0){
    G4cout.precision(5);       
    G4cout << " Projected Range= " << G4BestUnit(fProjRange, "Length")
           << "   rms= "             << G4BestUnit(rms, "Length")
           << "\n" << G4endl;
  };

  G4double ekin[100], dedxp[100], dedxmp[100], tdedxp[100], tdedxmp[100], 
    xsp[100], xsmp[100];
  G4EmCalculator calc;
  //calc.SetVerbose(2);
  G4int i;
  for(i = 0; i < 100; ++i) {
    ekin[i]   = std::pow(10., 0.1*G4double(i)) * keV;
    dedxp[i]  = calc.GetDEDX(ekin[i], proton, material);
    xsp[i]    = calc.GetCrossSectionPerVolume(ekin[i], proton, "hIoni", 
                                              material);
    tdedxp[i] = calc.ComputeElectronicDEDX(ekin[i], proton, material);
    dedxmp[i] = calc.GetDEDX(ekin[i], part, material);
    xsmp[i]   = calc.GetCrossSectionPerVolume(ekin[i], part, "mplIoni", 
                                              material);
    tdedxmp[i] = calc.ComputeElectronicDEDX(ekin[i], part, material);
  }

  if(GetVerbose() > 0){
    G4int prec = G4cout.precision(3);
    G4cout<<"##################################################################"
          << G4endl;
    G4cout<< "### Stopping Powers and Cross Sections" << G4endl;
    G4cout<<"##################################################################"
          << G4endl;
               
    G4cout<<"# N   E(MeV)  p_dEdx(MeV/mm) mpl_dEdx(MeV/mm)           xs(1/mm)" 
          << G4endl;
    G4cout<<"               restr    tot      restr    tot          p      mpl" 
          << G4endl;
    G4cout<<"##################################################################"
          << G4endl;
    for(i=0; i<100; ++i) {
      G4cout << std::setw(2) << i << "." << std::setw(9) << ekin[i] 
             << std::setw(8) << dedxp[i] 
             << std::setw(8) << tdedxp[i]
             << std::setw(9) << dedxmp[i] 
             << std::setw(9) << tdedxmp[i]
             << std::setw(10) << xsp[i]
             << std::setw(10) << xsmp[i]
             << G4endl;
    }
    G4cout.precision(prec);
    G4cout<<"##################################################################"
          << G4endl;
  }

  // normalize histogram
  G4double fac = (mm/MeV) / (nEvents * binLength);
  fAnalysisManager->ScaleH1(1,fac);

  for(i=0; i<100; ++i) {
    G4double e = std::log10(ekin[i] / MeV) + 0.05;
    fAnalysisManager->FillH1(2, e, tdedxp[i]);
    fAnalysisManager->FillH1(3, e, tdedxmp[i]);
    fAnalysisManager->FillH1(4, e, 
                   std::log10(calc.GetRange(ekin[i],"proton",matName)/mm));
    fAnalysisManager->FillH1(5, e, 
                   std::log10(calc.GetRange(ekin[i],"monopole",matName)/mm));
    fAnalysisManager->FillH1(6, e, dedxp[i]);
    fAnalysisManager->FillH1(7, e, dedxmp[i]);
    fAnalysisManager->FillH1(8, e, xsp[i]);
    fAnalysisManager->FillH1(9, e, xsmp[i]);
  }
}   

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void Run::FillHisto(G4int histoId, G4double v1, G4double v2)
{
  fAnalysisManager->FillH1(histoId, v1, v2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
