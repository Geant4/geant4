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
/// \file medical/electronScattering/src/Run.cc
/// \brief Implementation of the Run class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
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
#include <stdexcept>      // std::out_of_range

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det)
:  G4Run(), 
  fDetector(det),fParticle(0),
  fEnergy(0)
{
  InitFluence();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy )
{
  fParticle = particle;
  fEnergy = energy;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
   const Run* localRun = static_cast<const Run*>(run);

   // pass information about primary particle
   fParticle = localRun->fParticle;
   fEnergy     = localRun->fEnergy;

   // In MT all threads have the same fNbBins and fDr value
   fNbBins = localRun->fNbBins ;
   fDr     = localRun->fDr;

   // Some value by value 
   
  for (unsigned  int i = 0; i < localRun->fluence.size(); ++i)
  {
    try { fluence[i]+=localRun->fluence[i]; }
    catch(const std::out_of_range&) 
      { fluence.push_back(localRun->fluence[i]); }
  }

  for (unsigned  int i = 0; i < localRun->fluence1.size(); ++i)
  {
    try { fluence1[i]+=localRun->fluence1[i]; }
    catch(const std::out_of_range&) 
      { fluence1.push_back(localRun->fluence1[i]); }
  }

  for (unsigned  int i = 0; i < localRun->fluence2.size(); ++i)
  {
    try { fluence2[i]+=localRun->fluence2[i]; }
    catch(const std::out_of_range&) 
      { fluence2.push_back(localRun->fluence2[i]); }
  }
    
  for (unsigned  int i = 0; i < localRun->fNbEntries.size(); ++i)
  {
    try { fNbEntries[i]+=localRun->fNbEntries[i]; }
    catch(const std::out_of_range&) 
      { fNbEntries.push_back(localRun->fNbEntries[i]); }
  }

  G4Run::Merge(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void Run::EndOfRun()
{
  

  if (numberOfEvent == 0) return;
      
  //Scatter foil
  
  G4Material* material = fDetector->GetMaterialScatter();
  G4double length  = fDetector->GetThicknessScatter();
  G4double density = material->GetDensity();
  G4String partName = fParticle->GetParticleName();

  G4cout << "\n ======================== run summary ======================\n";

  G4int prec = G4cout.precision(3);
  
  G4cout << "\n The run was " << numberOfEvent << " " << partName << " of "
         << G4BestUnit(fEnergy,"Energy") << " through " 
         << G4BestUnit(length,"Length") << " of "
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;

  G4cout.precision(prec);

  // normalize histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
  G4double fac = 1./(double(numberOfEvent));
  analysisManager->ScaleH1(1,fac);
  analysisManager->ScaleH1(2,fac);
  analysisManager->ScaleH1(3,fac);
  analysisManager->ScaleH1(5,fac);
  analysisManager->ScaleH1(6,fac);
  
  ComputeFluenceError();
  PrintFluence(numberOfEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::InitFluence()
{
  // create Fluence histo in any case
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4int ih = 4;
  analysisManager->SetH1(ih, 120, 0*mm, 240*mm, "mm");
    
  //construct vectors for fluence distribution
  
  fNbBins = analysisManager->GetH1Nbins(ih);
  fDr = analysisManager->GetH1Width(ih);
  fluence.resize(fNbBins, 0.); 
  fluence1.resize(fNbBins, 0.);   
  fluence2.resize(fNbBins, 0.);
  fNbEntries.resize(fNbBins, 0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SumFluence(G4double r, G4double fl)
{
  G4int ibin = (int)(r/fDr);
  if (ibin >= fNbBins) return;
  fNbEntries[ibin]++;
  fluence[ibin]  += fl;
  fluence2[ibin] += fl*fl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::ComputeFluenceError()
{
  //compute rms
 
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
  
  G4double rnorm(4*mm), radius(0.), fnorm(0.), fnorm2(0.);
  G4int inorm = -1;
  do {
   inorm++; radius += fDr; fnorm += fluence[inorm]; fnorm2 += fluence2[inorm];
  } while (radius < rnorm);
  fnorm  /= (inorm+1);
  fnorm2 /= (inorm+1);  
    
  G4double ratio, error;
  G4double scale = 1./fnorm;
  G4double err0 = fnorm2/fnorm, err1 = 0.;
  
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
     G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();     
     analysisManager->FillH1(4,rmean,ratio);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <fstream>

void Run::PrintFluence(G4int TotEvents)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();   
  G4String name = analysisManager->GetFileName(); 
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

