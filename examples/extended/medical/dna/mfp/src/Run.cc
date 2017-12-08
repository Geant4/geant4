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

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(const DetectorConstruction* detector)
: G4Run(),
  fDetector(detector),
  fParticle(0), fEkin(0.),  
  fTotalCount(0), fSumTrack(0.), fSumTrack2(0.), fEnTransfer(0.)
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

  // map: processes count
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
  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  G4int prec = G4cout.precision(2);
  
  // run conditions  
  G4Material* material = fDetector->GetAbsorMaterial();
  G4double density  = material->GetDensity();       
  G4String partName = fParticle->GetParticleName();
  
  G4cout << 
   "\n ======================== run summary =====================\n";  
  G4cout 
    << "\n The run is " << numberOfEvent << " "<< partName << " of "
    << G4BestUnit(fEkin,"Energy") << " through a sphere of radius "
    << G4BestUnit(fDetector->GetAbsorRadius(),"Length") << "of "
    << material->GetName() << " (density: " 
    << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;    

  if (numberOfEvent == 0) {
    G4cout.setf(mode,std::ios::floatfield);
    G4cout.precision(prec);  
    return;
  }
                    
  // frequency of processes
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
           << "a radius of "
           << G4BestUnit(fDetector->GetAbsorRadius(),"Length") << " of "
           << material->GetName() << " : " << survive << G4endl;
  }

  if (fTotalCount == 0) fTotalCount = 1;   //force printing anyway

  // compute mean free path and related quantities
  G4double MeanFreePath = fSumTrack /fTotalCount;     
  G4double MeanTrack2   = fSumTrack2/fTotalCount;     
  G4double rmsBis = 
    std::sqrt(std::fabs(MeanTrack2 - MeanFreePath*MeanFreePath));
  G4double CrossSection = 1./MeanFreePath;     
  G4double massicMFP = MeanFreePath*density;
  G4double massicCS  = 1./massicMFP;

  G4cout << "\n\n MeanFreePath:\t"   << G4BestUnit(MeanFreePath,"Length")
         << " +- "                   << G4BestUnit(rmsBis,"Length")
         << "\t\t\tmassic: "         << G4BestUnit(massicMFP, "Mass/Surface")
         << "\n CrossSection:\t"     << CrossSection*cm << " cm^-1 "
         << "\t\t\tmassic: "         << G4BestUnit(massicCS, "Surface/Mass")
         << G4endl;
         
  // compute energy transfer coefficient
  G4double MeanTransfer   = fEnTransfer/fTotalCount;
  G4double massTransfCoef = massicCS*MeanTransfer/fEkin;
   
  G4cout << "\n mean energy of charged secondaries: " 
         << G4BestUnit(MeanTransfer, "Energy")
         << "\tmass_energy_transfer coef: "          
         << G4BestUnit(massTransfCoef, "Surface/Mass")
         << G4endl;       

  //output file
  //
  FILE* myFile;
  myFile=fopen("mfp.txt","a");
  fprintf(myFile,"%e %e %e \n",
     fEkin/eV,
     MeanFreePath/nm,
     rmsBis/nm);
  fclose(myFile);

  // remove all contents in fProcCounter 
  fProcCounter.clear();

  //reset default formats
  //
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);

}
