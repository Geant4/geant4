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
/// \file medical/electronScattering/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin,
                     HistoManager* histo)
:fDetector(det), fPrimary(kin), fHistoManager(histo)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  fHistoManager->book();
  
  InitFluence();  

  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int TotNbofEvents = aRun->GetNumberOfEvent();
  if (TotNbofEvents == 0) return;
      
  //Scatter foil
  //
  G4Material* material = fDetector->GetMaterialScatter();
  G4double length  = fDetector->GetThicknessScatter();
  G4double density = material->GetDensity();
   
  G4ParticleDefinition* particle = fPrimary->GetParticleGun()
                                          ->GetParticleDefinition();
  G4String partName = particle->GetParticleName();
  G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();

  G4cout << "\n ======================== run summary ======================\n";

  G4int prec = G4cout.precision(3);
  
  G4cout << "\n The run was " << TotNbofEvents << " " << partName << " of "
         << G4BestUnit(energy,"Energy") << " through " 
         << G4BestUnit(length,"Length") << " of "
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;

  G4cout.precision(prec);

  // normalize histograms
  //
  G4double fac = 1./(double(TotNbofEvents));
  fHistoManager->Normalize(1,fac);
  fHistoManager->Normalize(2,fac);
  fHistoManager->Normalize(3,fac);
  
  ComputeFluenceError();
  PrintFluence(TotNbofEvents);

  // save histograms
  fHistoManager->save();

  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::InitFluence()
{
  // create Fluence histo in any case
  //
  G4int ih = 4;
  if (!fHistoManager->HistoExist(ih))
    fHistoManager->SetHisto(ih, 120, 0*mm, 240*mm, "mm");
    
  //construct vectors for fluence distribution
  //
  fNbBins = fHistoManager->GetNbins(ih);
  fDr = fHistoManager->GetBinWidth(ih);
  fluence.resize(fNbBins, 0.); 
  fluence1.resize(fNbBins, 0.);   
  fluence2.resize(fNbBins, 0.);
  fNbEntries.resize(fNbBins, 0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SumFluence(G4double r, G4double fl)
{
  G4int ibin = (int)(r/fDr);
  if (ibin >= fNbBins) return;
  fNbEntries[ibin]++;
  fluence[ibin]  += fl;
  fluence2[ibin] += fl*fl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::ComputeFluenceError()
{
  //compute rms
  //
  G4double ds,variance,rms;
  G4double rmean = -0.5*fDr;

  for (G4int bin=0; bin<fNbBins; bin++) {
     rmean += fDr;  
     ds = twopi*rmean*fDr;
     fluence[bin] /= ds;
     fluence2[bin] /= (ds*ds);
     variance = 0.;
     if (fNbEntries[bin] > 0)     
       variance = fluence2[bin] - (fluence[bin]*fluence[bin])/fNbEntries[bin];
     rms = 0.;
     if(variance > 0.) rms = std::sqrt(variance);
     fluence2[bin] = rms;
  }
   
  //normalize to first bins, compute error and fill histo
  //
  G4double rnorm(4*mm), radius(0.), fnorm(0.), fnorm2(0.);
  G4int inorm = -1;
  do {
   inorm++; radius += fDr; fnorm += fluence[inorm]; fnorm2 += fluence2[inorm];
  } while (radius < rnorm);
  fnorm  /= (inorm+1);
  fnorm2 /= (inorm+1);  
  //  
  G4double ratio, error;
  G4double scale = 1./fnorm;
  G4double err0 = fnorm2/fnorm, err1 = 0.;
  //
  rmean = -0.5*fDr;
    
  for (G4int bin=0; bin<fNbBins; bin++) {
     ratio = fluence[bin]*scale;
     error = 0.;
     if (ratio > 0.) {
       err1  = fluence2[bin]/fluence[bin]; 
       error = ratio*std::sqrt(err1*err1 + err0*err0);
     }
     fluence1[bin] = ratio;
     fluence2[bin] = error;
     rmean += fDr;
     fHistoManager->FillHisto(4,rmean,ratio);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <fstream>

void RunAction::PrintFluence(G4int TotEvents)
{
  G4String name = fHistoManager->GetFileName(); 
  G4String fileName = name + ".ascii";
  std::ofstream File(fileName, std::ios::out);
  
  std::ios::fmtflags mode = File.flags();  
  File.setf( std::ios::scientific, std::ios::floatfield );
  G4int prec = File.precision(3);
      
  File << "  Fluence density distribution \n " 
       << "\n  ibin \t radius (mm) \t Nb \t fluence\t norma fl\t rms/nfl (%) \n"
       << G4endl;

  G4double rmean = -0.5*fDr;    
  for (G4int bin=0; bin<fNbBins; bin++) {
     rmean +=fDr;
     G4double error = 0.;
     if (fluence1[bin] > 0.) error =  100*fluence2[bin]/fluence1[bin];
     File << "  " << bin << "\t " << rmean/mm << "\t " << fNbEntries[bin]
          << "\t " << fluence[bin]/double(TotEvents) << "\t " << fluence1[bin] 
          << "\t " << error
          << G4endl;          
  }
    
  // restaure default formats
  File.setf(mode,std::ios::floatfield);
  File.precision(prec);         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

