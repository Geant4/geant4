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
/// \file electromagnetic/TestEm11/src/Run.cc
/// \brief Implementation of the Run class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "DetectorConstruction.hh"

#include "EventAction.hh"
#include "HistoManager.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Material.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* detector)
: G4Run(),
  fDetector(detector),
  fParticle(0), fEkin(0.),  
  fTrackLen(0.),  fTrackLen2(0.),
  fProjRange(0.), fProjRange2(0.),
  fNbOfSteps(0), fNbOfSteps2(0),
  fStepSize(0.),  fStepSize2(0.)
{
  for (G4int i=0; i<3; ++i) { fStatus[i] = 0; fTotEdep[i] = 0.; }
  fTotEdep[1] = joule;
  for (G4int i=0; i<kMaxAbsor; ++i) {
    fEdeposit[i] = 0.; fEmin[i] = joule; fEmax[i] = 0.;
    fCsdaRange[i] = 0.; fXfrontNorm[i] = 0.;
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimary (G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin     = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddEdep (G4int i, G4double e)        
{
  if (e > 0.) {
    fEdeposit[i]  += e;
    if (e < fEmin[i]) fEmin[i] = e;
    if (e > fEmax[i]) fEmax[i] = e;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddTotEdep (G4double e)        
{
  if (e > 0.) {
    fTotEdep[0]  += e;
    if (e < fTotEdep[1]) fTotEdep[1] = e;
    if (e > fTotEdep[2]) fTotEdep[2] = e;
  }
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
      
void Run::AddTrackStatus (G4int i)    
{
  fStatus[i]++ ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
void Run::SetCsdaRange (G4int i, G4double value) 
{
  fCsdaRange[i] = value; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
void Run::SetXfrontNorm (G4int i, G4double value) 
{
  fXfrontNorm[i] = value; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
                                      
G4double Run::GetCsdaRange (G4int i) 
{
  return fCsdaRange[i];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
G4double Run::GetXfrontNorm (G4int i) 
{
  return fXfrontNorm[i];
}
       
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);
  
  // pass information about primary particle
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;

  // accumulate sums
  fTrackLen   += localRun->fTrackLen;  
  fTrackLen2  += localRun->fTrackLen2;
  fProjRange  += localRun->fProjRange; 
  fProjRange2 += localRun->fProjRange2;
  fNbOfSteps  += localRun->fNbOfSteps ;
  fNbOfSteps2 += localRun->fNbOfSteps2;
  fStepSize   += localRun->fStepSize;  
  fStepSize2  += localRun->fStepSize2;
  
  G4int nbOfAbsor = fDetector->GetNbOfAbsor();
  for (G4int i=1; i<=nbOfAbsor; ++i) {
    fEdeposit[i]  += localRun->fEdeposit[i];
    fCsdaRange[i]  = localRun->fCsdaRange[i]; 
    fXfrontNorm[i] = localRun->fXfrontNorm[i];
    // min, max
    G4double min,max;
    min = localRun->fEmin[i]; max = localRun->fEmax[i];
    if (fEmin[i] > min) fEmin[i] = min;
    if (fEmax[i] < max) fEmax[i] = max;
  } 
   
  for (G4int i=0; i<3; ++i)  fStatus[i] += localRun->fStatus[i];
  
  // total Edep
  fTotEdep[0] += localRun->fTotEdep[0];
  G4double min,max;
  min = localRun->fTotEdep[1]; max = localRun->fTotEdep[2];
if (fTotEdep[1] > min) fTotEdep[1] = min;
if (fTotEdep[2] < max) fTotEdep[2] = max;

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
  G4int nbOfAbsor   = fDetector->GetNbOfAbsor();
  
  G4cout << "\n ======================== run summary =====================\n";
  G4cout 
    << "\n The run is " << numberOfEvent << " "<< partName << " of "
    << G4BestUnit(fEkin,"Energy") 
    << " through "  << nbOfAbsor << " absorbers: \n";
  for (G4int i=1; i<= nbOfAbsor; i++) {
     G4Material* material = fDetector->GetAbsorMaterial(i);
     G4double thickness = fDetector->GetAbsorThickness(i);
     G4double density = material->GetDensity();
     G4cout << std::setw(5) << i
            << std::setw(10) << G4BestUnit(thickness,"Length") << " of "
            << material->GetName() << " (density: " 
            << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
  }         

  if (numberOfEvent == 0) {
    G4cout.setf(mode,std::ios::floatfield);
    G4cout.precision(prec);  
    return;
  }
  
  G4cout.precision(3);
  G4double rms (0);
  
  for (G4int i=1; i<= nbOfAbsor; i++) {
     fEdeposit[i] /= numberOfEvent;

     G4cout 
       << "\n Edep in absorber " << i << " = " 
       << G4BestUnit(fEdeposit[i],"Energy")
       << "\t(" << G4BestUnit(fEmin[i], "Energy")
       << "-->" << G4BestUnit(fEmax[i], "Energy")
       << ")";
  }
  G4cout << G4endl;

  if (nbOfAbsor > 1) {
    fTotEdep[0] /= numberOfEvent;
    G4cout 
      << "\n Edep in all absorbers = " << G4BestUnit(fTotEdep[0],"Energy")
      << "\t(" << G4BestUnit(fTotEdep[1], "Energy")
      << "-->" << G4BestUnit(fTotEdep[2], "Energy")
      << ")" << G4endl;
  }
  
  //compute track length of primary track
  //
  fTrackLen /= numberOfEvent; fTrackLen2 /= numberOfEvent;
  rms = fTrackLen2 - fTrackLen*fTrackLen;        
  if (rms>0.) rms = std::sqrt(rms); else rms = 0.;

  G4cout.precision(3);       
  G4cout 
    << "\n Track length of primary track = " << G4BestUnit(fTrackLen,"Length")
    << " +- "                                << G4BestUnit( rms,"Length");
    
  //compare with csda range
  //
  G4int NbOfAbsor = fDetector->GetNbOfAbsor();
  if (NbOfAbsor == 1) {
    G4cout 
     << "\n Range from EmCalculator = " << G4BestUnit(fCsdaRange[1],"Length")
     << " (from full dE/dx)" << G4endl;
  }
                     
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
  G4cout << "\n Nb of steps of primary track  = " << fNbSteps << " +- " << rms;
    
  fStepSize /= numberOfEvent; fStepSize2 /= numberOfEvent;
  rms = fStepSize2 - fStepSize*fStepSize;        
  if (rms>0.) rms = std::sqrt(rms); else rms = 0.;

  G4cout.precision(3);       
  G4cout 
    << "\t Step size= " << G4BestUnit(fStepSize,"Length")
    << " +- "           << G4BestUnit( rms,"Length")
    << G4endl;
    
  //transmission coefficients
  //
  G4double absorbed  = 100.*fStatus[0]/dNofEvents;
  G4double transmit  = 100.*fStatus[1]/dNofEvents;
  G4double reflected = 100.*fStatus[2]/dNofEvents;  

  G4cout.precision(2);       
  G4cout 
    << "\n absorbed = "  << absorbed  << " %"
    << "   transmit = "  << transmit  << " %"
    << "   reflected = " << reflected << " %" << G4endl;

    // normalize histograms of longitudinal energy profile
    //
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4int ih = 1;
    G4double binWidth = analysisManager->GetH1Width(ih)
                       *analysisManager->GetH1Unit(ih);
    G4double fac = (1./(numberOfEvent*binWidth))*(mm/MeV);
    analysisManager->ScaleH1(ih,fac);
    
    ih = 8;
    binWidth = analysisManager->GetH1Width(ih);
    fac = (1./(numberOfEvent*binWidth))*(g/(MeV*cm2));
    analysisManager->ScaleH1(ih,fac);
    
   // reset default formats
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
