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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 37 (2010) 4692-4708
// Phys. Med. 31 (2015) 861-874
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file medical/dna/svalue/src/Run.cc
/// \brief Implementation of the Run class

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "HistoManager.hh"
#include "MyFile.hh"

#ifdef MYFILE
 #include "MyPrimaryGeneratorActionFromFile.hh"
#else
 #include "PrimaryGeneratorAction.hh"
#endif

#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(const DetectorConstruction* detector)
: G4Run(),
  fDetector(detector),
  fParticle(0), fEkin(0.),  
  fCytoEdeposit(0.),  fCytoEdeposit2(0.),
  fNuclEdeposit(0.),  fNuclEdeposit2(0.),
  fTrackLen(0.),  fTrackLen2(0.),
  fProjRange(0.), fProjRange2(0.),
  fNbOfSteps(0), fNbOfSteps2(0),
  fStepSize(0.),  fStepSize2(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimary (G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin     = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddCytoEdep (G4double e)        
{
  fCytoEdeposit  += e;
  fCytoEdeposit2 += e*e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddNuclEdep (G4double e)        
{
  fNuclEdeposit  += e;
  fNuclEdeposit2 += e*e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
void Run::AddTrackLength (G4double t) 
{
  fTrackLen  += t;
  fTrackLen2 += t*t;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
void Run::AddProjRange (G4double x) 
{
  fProjRange  += x;
  fProjRange2 += x*x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
void Run::AddStepSize (G4int nb, G4double st)
{
  fNbOfSteps  += nb; 
  fNbOfSteps2 += nb*nb;
  fStepSize   += st ; 
  fStepSize2  += st*st;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);
  
  // pass information about primary particle
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;

  // accumulate sums
  
  fCytoEdeposit   += localRun->fCytoEdeposit;
  fCytoEdeposit2  += localRun->fCytoEdeposit2;
  fNuclEdeposit   += localRun->fNuclEdeposit;
  fNuclEdeposit2  += localRun->fNuclEdeposit2;
  
  fTrackLen   += localRun->fTrackLen;  
  fTrackLen2  += localRun->fTrackLen2;
  fProjRange  += localRun->fProjRange; 
  fProjRange2 += localRun->fProjRange2;
  fNbOfSteps  += localRun->fNbOfSteps ;
  fNbOfSteps2 += localRun->fNbOfSteps2;
  fStepSize   += localRun->fStepSize;  
  fStepSize2  += localRun->fStepSize2;

  G4Run::Merge(run); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun() 
{
  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  G4int prec = G4cout.precision(2);
  
  //run conditions  
  //
  G4String partName = fParticle->GetParticleName();
  
  G4cout << "\n ======================== run summary =====================\n";  
  G4cout 
    << "\n The run is " << numberOfEvent << " "<< partName << " of "
    << G4BestUnit(fEkin,"Energy") << G4endl;    

  if (numberOfEvent == 0) {
    G4cout.setf(mode,std::ios::floatfield);
    G4cout.precision(prec);  
    return;
  }
      
  //
  
  fCytoEdeposit /= numberOfEvent; fCytoEdeposit2 /= numberOfEvent;
  G4double rmsCyto = fCytoEdeposit2 - fCytoEdeposit*fCytoEdeposit;        
  if (rmsCyto>0.) rmsCyto = std::sqrt(rmsCyto); else rmsCyto = 0.;

  G4cout.precision(3);       
  G4cout 
    << "\n Total Energy deposited in cytoplasm = " << G4BestUnit(fCytoEdeposit,"Energy")
    << " +- "                                << G4BestUnit( rmsCyto,"Energy")
    << G4endl;
                    
  G4double sValueCyto=fCytoEdeposit/fDetector->GetCytoMass();
  G4double rmsSValueCyto=rmsCyto/fDetector->GetCytoMass();
  
  G4cout.precision(3);       
  G4cout 
    << "\n S value for cytoplasm (C<-C) = " << sValueCyto/gray << " Gy/Bq.s "
    << " +- "                                << rmsSValueCyto/gray 
    <<  " Gy/Bq.s "
    << G4endl;
              
  //
  
  fNuclEdeposit /= numberOfEvent; fNuclEdeposit2 /= numberOfEvent;
  G4double rmsNucl = fNuclEdeposit2 - fNuclEdeposit*fNuclEdeposit;        
  if (rmsNucl>0.) rmsNucl = std::sqrt(rmsNucl); else rmsNucl = 0.;

  G4cout.precision(3);       
  G4cout 
    << "\n Total Energy deposited in nucleus = " << G4BestUnit(fNuclEdeposit,"Energy")
    << " +- "                                << G4BestUnit( rmsNucl,"Energy")
    << G4endl;
                    
  G4double sValueNucl=fNuclEdeposit/fDetector->GetNuclMass();
  G4double rmsSValueNucl=rmsNucl/fDetector->GetNuclMass();
  
  G4cout.precision(3);       
  G4cout 
    << "\n S value for nucleus (N<-C) = " << sValueNucl/gray << " Gy/Bq.s "
    << " +- "                                << rmsSValueNucl/gray 
    <<  " Gy/Bq.s "
    << G4endl;

  //compute track length of primary track
  //
  fTrackLen /= numberOfEvent; fTrackLen2 /= numberOfEvent;
  G4double rms = fTrackLen2 - fTrackLen*fTrackLen;        
  if (rms>0.) rms = std::sqrt(rms); else rms = 0.;

  G4cout.precision(3);       
  G4cout 
    << "\n Track length of primary track = " << G4BestUnit(fTrackLen,"Length")
    << " +- "                                << G4BestUnit( rms,"Length");
    
  //compute projected range of primary track
  //
  fProjRange /= numberOfEvent; fProjRange2 /= numberOfEvent;
  rms = fProjRange2 - fProjRange*fProjRange;        
  if (rms>0.) rms = std::sqrt(rms); else rms = 0.;
   
  G4cout 
    << "\n Projected range               = " << G4BestUnit(fProjRange,"Length")
    << " +- "                                << G4BestUnit( rms,"Length")    
    << G4endl;
    
  //nb of steps and step size of primary track
  //
  G4double dNofEvents = double(numberOfEvent);
  G4double fNbSteps  = fNbOfSteps/dNofEvents, 
           fNbSteps2 = fNbOfSteps2/dNofEvents;
  rms = fNbSteps2 - fNbSteps*fNbSteps;       
  if (rms>0.) rms = std::sqrt(rms); else rms = 0.;

  G4cout.precision(2);       
  G4cout << "\n Nb of steps of primary track  = " << fNbSteps << " +- " << rms
  << G4endl;
    
  fStepSize /= numberOfEvent; fStepSize2 /= numberOfEvent;
  rms = fStepSize2 - fStepSize*fStepSize;        
  if (rms>0.) rms = std::sqrt(rms); else rms = 0.;

  G4cout.precision(3);       
  G4cout 
    << "\n Step size                     = " << G4BestUnit(fStepSize,"Length")
    << " +- "           << G4BestUnit( rms,"Length")
    << G4endl;

  // normalize histograms of longitudinal energy profile
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4int ih = 1;
  G4double binWidth = analysisManager->GetH1Width(ih);
  G4double fac = (1./(numberOfEvent*binWidth))*(mm/MeV);
  analysisManager->ScaleH1(ih,fac);
    
  // reset default formats
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);
  
  //output file
  FILE *myFile;
  myFile = fopen ("s.txt","a");
  fprintf (myFile, "%e %e %e %e %e %e %e \n", 
    fDetector->GetNuclRadius()/nm,
    fDetector->GetCytoThickness()/nm,
    fEkin/eV,
    sValueCyto/gray,
    rmsSValueCyto/gray,
    sValueNucl/gray,
    rmsSValueNucl/gray
    );
  fclose (myFile);
  
}
