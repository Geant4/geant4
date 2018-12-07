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
/// \file electromagnetic/TestEm13/src/Run.cc
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
  fParticle(0), fEkin(0.)
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
  G4double tickness = fDetector->GetSize();
     
  G4cout << "\n ======================== run summary ======================\n";
  G4cout << "\n The run is: " << numberOfEvent << " " << partName << " of "
         << G4BestUnit(fEkin,"Energy") << " through " 
         << G4BestUnit(tickness,"Length") << " of "
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;

  //frequency of processes
  G4int totalCount = 0;
  G4int survive = 0;  
  G4cout << "\n Process calls frequency --->";
  std::map<G4String,G4int>::iterator it;  
  for (it = fProcCounter.begin(); it != fProcCounter.end(); it++) {
     G4String procName = it->first;
     G4int    count    = it->second;
     totalCount += count; 
     G4cout << "\t" << procName << " = " << count;
     if (procName == "Transportation") survive = count;
  }
  G4cout << G4endl;

  if (totalCount == 0) { G4cout.precision(dfprec);   return;};  
  G4double ratio = double(survive)/totalCount;

  G4cout << "\n Nb of incident particles unaltered after "
         << G4BestUnit(tickness,"Length") << " of "
         << material->GetName() << " : " << survive 
         << " over " << totalCount << " incident particles."
         << "  Ratio = " << 100*ratio << " %" << G4endl;
  
  if (ratio == 0.) return;
  
  //compute cross section and related quantities
  //
  G4double CrossSection = - std::log(ratio)/tickness;     
  G4double massicCS  = CrossSection/density;
   
  G4cout << " ---> CrossSection per volume:\t" << CrossSection*cm << " cm^-1 "
         << "\tCrossSection per mass: " << G4BestUnit(massicCS, "Surface/Mass")
         << G4endl;

  //check cross section from G4EmCalculator
  //
  G4cout << "\n Verification from G4EmCalculator: \n"; 
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
    if (procName != "Transportation")
      G4cout << "\t" << procName << "= " 
             << G4BestUnit(massSigma, "Surface/Mass");
  }             
  G4cout << "\ttotal= " 
         << G4BestUnit(sumc, "Surface/Mass") << G4endl;
         
  //expected ratio of transmitted particles
  G4double Ratio = std::exp(-sumc*density*tickness);
  G4cout << "\tExpected ratio of transmitted particles= " 
         << 100*Ratio << " %" << G4endl;         
                        
  // remove all contents in fProcCounter 
  fProcCounter.clear();
  
  //restore default format
  G4cout.precision(dfprec);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
