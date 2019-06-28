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
/// \file electromagnetic/TestEm1/src/Run.cc
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

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(const DetectorConstruction* det)
: G4Run(),
  fDetector(det), 
  fParticle(0), fEkin(0.),
  fNbOfTraks0(0), fNbOfTraks1(0),
  fNbOfSteps0(0), fNbOfSteps1(0),
  fEdep(0.), fNIEL(0.),
  fTrueRange(0.), fTrueRange2(0.),
  fProjRange(0.), fProjRange2(0.),  
  fTransvDev(0.), fTransvDev2(0.)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimary(const G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin = energy;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::CountProcesses(const G4String& procName) 
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

  // accumulate sums
  //
  fNbOfTraks0 += localRun->fNbOfTraks0;  
  fNbOfTraks1 += localRun->fNbOfTraks1;  
  fNbOfSteps0 += localRun->fNbOfSteps0;
  fNbOfSteps1 += localRun->fNbOfSteps1;   
  fEdep       += localRun->fEdep;  
  fNIEL       += localRun->fNIEL;  
  fTrueRange  += localRun->fTrueRange;
  fTrueRange2 += localRun->fTrueRange2;
  fProjRange  += localRun->fProjRange;
  fProjRange2 += localRun->fProjRange2;
  fTransvDev  += localRun->fTransvDev;
  fTransvDev2 += localRun->fTransvDev2;  
      
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
  G4int prec = 5, wid = prec + 2;  
  G4int dfprec = G4cout.precision(prec);
  
  //run condition
  //        
  G4String partName    = fParticle->GetParticleName();    
  const G4Material* material = fDetector->GetMaterial();
  G4double density     = material->GetDensity();
     
  G4cout << "\n ======================== run summary ======================\n";
  G4cout << "\n The run is: " << numberOfEvent << " " << partName << " of "
         << G4BestUnit(fEkin,"Energy") << " through " 
         << G4BestUnit(fDetector->GetSize(),"Length") << " of "
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;

  if (numberOfEvent == 0) { G4cout.precision(dfprec);   return;}   

  G4double dNbOfEvents = (G4double)numberOfEvent;  
  G4cout << "\n Total energy deposit:   " 
         << G4BestUnit(fEdep/dNbOfEvents, "Energy") << G4endl;
  G4cout << " NIEL energy calculated: " 
         << G4BestUnit(fNIEL/dNbOfEvents, "Energy") << G4endl;
               
  //nb of tracks and steps per event
  //           
  G4cout << "\n Nb tracks/event"
         << "   neutral: " << std::setw(wid) << fNbOfTraks0/dNbOfEvents
         << "   charged: " << std::setw(wid) << fNbOfTraks1/dNbOfEvents
         << "\n Nb  steps/event"
         << "   neutral: " << std::setw(wid) << fNbOfSteps0/dNbOfEvents
         << "   charged: " << std::setw(wid) << fNbOfSteps1/dNbOfEvents
         << G4endl;
        
  //frequency of processes
  //
  G4cout << "\n Process calls frequency :" << G4endl;
  G4int index = 0;  
  std::map<G4String,G4int>::iterator it;         
  for (it = fProcCounter.begin(); it != fProcCounter.end(); it++) {
     G4String procName = it->first;
     G4int    count    = it->second;
     G4String space = " "; if (++index%3 == 0) space = "\n";
     G4cout << " " << std::setw(20) << procName << "="<< std::setw(7) << count
            << space;
  }     
  G4cout << G4endl;
        
  //compute true and projected ranges, and transverse dispersion
  //
  fTrueRange /= numberOfEvent; fTrueRange2 /= numberOfEvent;
  G4double trueRms = fTrueRange2 - fTrueRange*fTrueRange;        
  if (trueRms>0.) trueRms = std::sqrt(trueRms); else trueRms = 0.;
        
  fProjRange /= numberOfEvent; fProjRange2 /= numberOfEvent;
  G4double projRms = fProjRange2 - fProjRange*fProjRange;        
  if (projRms>0.) projRms = std::sqrt(projRms); else projRms = 0.;
         
  fTransvDev /= 2*numberOfEvent; fTransvDev2 /= 2*numberOfEvent;
  G4double trvsRms = fTransvDev2 - fTransvDev*fTransvDev;        
  if (trvsRms>0.) trvsRms = std::sqrt(trvsRms); else trvsRms = 0.;
   
  //compare true range with csda range from PhysicsTables
  //  
  G4EmCalculator emCalculator;
  G4double rangeTable = 0.;
  if (fParticle->GetPDGCharge() != 0.)
    rangeTable = emCalculator.GetCSDARange(fEkin,fParticle,material);
        
  G4cout << "\n---------------------------------------------------------\n";
  G4cout << " Primary particle : " ;
  G4cout << "\n true Range = " << G4BestUnit(fTrueRange,"Length")
         << "   rms = "        << G4BestUnit(trueRms,  "Length");

  G4cout << "\n proj Range = " << G4BestUnit(fProjRange,"Length")
         << "   rms = "        << G4BestUnit(projRms,  "Length");
               
  G4cout << "\n proj/true  = " << fProjRange/fTrueRange;
                     
  G4cout << "\n transverse dispersion at end = " 
         << G4BestUnit(trvsRms,"Length");
          
  G4cout << "\n      mass true Range from simulation = " 
         << G4BestUnit(fTrueRange*density, "Mass/Surface")
         << "\n       from PhysicsTable (csda range) = " 
         << G4BestUnit(rangeTable*density, "Mass/Surface");        
  G4cout << "\n---------------------------------------------------------\n";
  G4cout << G4endl;
           
  // remove all contents in fProcCounter 
  fProcCounter.clear();
  
  //restore default format         
  G4cout.precision(dfprec);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
