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
// Phys. Med. 31 (2015) 861-874
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $Id$
//
/// \filemedical/dna/slowing/src/Run.cc
/// \brief Implementation of the Run class

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "HistoManager.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(const DetectorConstruction* detector)
: G4Run(),
  fDetector(detector),
  fParticle(0), fEkin(0.),  
  fEdeposit(0.),  fEdeposit2(0.)
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

void Run::AddEdep (G4double e)        
{
  fEdeposit  += e;
  fEdeposit2 += e*e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);
  
  // pass information about primary particle
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;

  // accumulate sums
  fEdeposit   += localRun->fEdeposit;
  fEdeposit2  += localRun->fEdeposit2;

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
  G4Material* material = fDetector->GetAbsorMaterial();
  G4double density  = material->GetDensity();       
    
  G4cout << "\n ======================== run summary =====================\n";  
  G4cout 
    << "\n The run is " << numberOfEvent << " "<< partName << " of "
    << G4BestUnit(fEkin,"Energy") << " through a volume of "
    << material->GetName() << " (density: " 
    << G4BestUnit(density,"Volumic Mass") << ") of mass " 
    << G4BestUnit(fDetector->GetAbsorMass(),"Mass") 
    << G4endl;    

  if (numberOfEvent == 0) {
    G4cout.setf(mode,std::ios::floatfield);
    G4cout.precision(prec);  
    return;
  }
      
  G4cout.precision(3);       
  G4cout 
    << "\n Total Energy deposited        = " << G4BestUnit(fEdeposit,"Energy")
    << G4endl;
                    
  /*
  G4double dose=fEdeposit/fDetector->GetAbsorMass();
  G4double rmsDose=rms/fDetector->GetAbsorMass();
  
  G4cout.precision(3);       
  G4cout 
    << "\n Dose                         = " << dose/gray << " Gy "
    << G4endl;
  */
              
  G4cout << G4endl;

  // normalize histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4double fac = 1./(numberOfEvent*(fEkin/eV));
  analysisManager->ScaleH1(1,fac);
  analysisManager->ScaleH1(2,fac);
  analysisManager->ScaleH1(3,fac);
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
