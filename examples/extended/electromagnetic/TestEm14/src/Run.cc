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
/// \file Run.cc
/// \brief Implementation of the Run class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmCalculator.hh"
#include "G4Gamma.hh"

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det)
: G4Run(),
  fDetector(det), 
  fParticle(0), fEkin(0.),
  fTotalCount(0), fSumTrack(0.), fSumTrack2(0.), fEnTransfer(0.)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin = energy;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::CountProcesses(G4String procName) 
{
  std::map<G4String,G4int>::iterator it = fProcCounter.find(procName);
  if ( it == fProcCounter.end()) {
    fProcCounter[procName] = 1;
  }
  else {
    fProcCounter[procName]++; 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SumTrack (G4double track)
{
  fTotalCount++;
  fSumTrack += track;
  fSumTrack2 += track*track;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SumeTransf (G4double energy)
{
  fEnTransfer += energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);

  // pass information about primary particle
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;
      
  //map: processes count
  std::map<G4String,G4int>::const_iterator it;
  for (it = localRun->fProcCounter.begin(); 
       it !=localRun->fProcCounter.end(); ++it) {
       
    G4String procName = it->first;
    G4int localCount  = it->second;
    if ( fProcCounter.find(procName) == fProcCounter.end()) {
      fProcCounter[procName] = localCount;
    }
    else {
      fProcCounter[procName] += localCount;
    }
  }
  
  fTotalCount += localRun->fTotalCount;
  fSumTrack   += localRun->fSumTrack;
  fSumTrack2  += localRun->fSumTrack2;
  fEnTransfer += localRun->fEnTransfer;
  
  G4Run::Merge(run); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun()
{
  G4int prec = 5;  
  G4int dfprec = G4cout.precision(prec);
  
  //run condition
  //        
  G4String partName    = fParticle->GetParticleName();    
  G4Material* material = fDetector->GetMaterial();
  G4double density     = material->GetDensity();
  G4double tickness    = fDetector->GetSize();
     
  G4cout << "\n ======================== run summary ======================\n";
  G4cout << "\n The run is: " << numberOfEvent << " " << partName << " of "
         << G4BestUnit(fEkin,"Energy") << " through " 
         << G4BestUnit(tickness,"Length") << " of "
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;

  //frequency of processes
  G4int survive = 0;  
  G4cout << "\n Process calls frequency --->";
  std::map<G4String,G4int>::iterator it;  
  for (it = fProcCounter.begin(); it != fProcCounter.end(); it++) {
     G4String procName = it->first;
     G4int    count    = it->second;
     G4cout << "\t" << procName << " = " << count;
     if (procName == "Transportation") survive = count;
  }

  if (survive > 0) {
    G4cout << "\n\n Nb of incident particles surviving after "
           << G4BestUnit(fDetector->GetSize(),"Length") << " of "
           << material->GetName() << " : " << survive << G4endl;
  }

  if (fTotalCount == 0) fTotalCount = 1;   //force printing anyway

  //compute mean free path and related quantities
  //
  G4double MeanFreePath = fSumTrack /fTotalCount;     
  G4double MeanTrack2   = fSumTrack2/fTotalCount;     
  G4double rms = std::sqrt(std::fabs(MeanTrack2 - MeanFreePath*MeanFreePath));
  G4double CrossSection = 1./MeanFreePath;     
  G4double massicMFP = MeanFreePath*density;
  G4double massicCS  = 1./massicMFP;

  G4cout << "\n\n MeanFreePath:\t"   << G4BestUnit(MeanFreePath,"Length")
         << " +- "                   << G4BestUnit( rms,"Length")
         << "\tmassic: "             << G4BestUnit(massicMFP, "Mass/Surface")
         << "\n CrossSection:\t"     << CrossSection*cm << " cm^-1 "
         << "\t\t\tmassic: "         << G4BestUnit(massicCS, "Surface/Mass")
         << G4endl;
         
  //compute energy transfer coefficient
  //
  G4double MeanTransfer   = fEnTransfer/fTotalCount;
  G4double massTransfCoef = massicCS*MeanTransfer/fEkin;
   
  G4cout << "\n mean energy of charged secondaries: " 
         << G4BestUnit(MeanTransfer, "Energy")
         << "\n     ---> mass_energy_transfer coef: "
         << G4BestUnit(massTransfCoef, "Surface/Mass")
         << G4endl;       

  //check cross section from G4EmCalculator
  //
  G4cout << "\n Verification : "
         << "crossSections from G4EmCalculator \n";

  G4EmCalculator emCalculator;
  G4double sumc = 0.0;  
  for (it = fProcCounter.begin(); it != fProcCounter.end(); it++) {
    G4String procName = it->first;      
    G4double massSigma = 
    emCalculator.GetCrossSectionPerVolume(fEkin,fParticle,
                                              procName,material)/density;
    if (fParticle == G4Gamma::Gamma())
       massSigma = 
       emCalculator.ComputeCrossSectionPerVolume(fEkin,fParticle,
                                              procName,material)/density;
    sumc += massSigma;
    G4cout << "   " << procName << "= " 
           << G4BestUnit(massSigma, "Surface/Mass");
  }
  G4cout << "   total= " 
         << G4BestUnit(sumc, "Surface/Mass") << G4endl;

  // remove all contents in fProcCounter 
  fProcCounter.clear();

  //restore default format
  G4cout.precision(dfprec);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
