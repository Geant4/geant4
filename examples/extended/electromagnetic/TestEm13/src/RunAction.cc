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
// $Id: RunAction.cc,v 1.8 2010-04-05 13:45:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"
#include "G4Gamma.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* prim)
  : detector(det), primary(prim)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
  
  G4int  prec = G4cout.precision(5);
    
  G4Material* material = detector->GetMaterial();
  G4double density  = material->GetDensity();
  G4double tickness = detector->GetSize();
   
  G4ParticleDefinition* particle = 
                            primary->GetParticleGun()->GetParticleDefinition();
  G4String Particle = particle->GetParticleName();    
  G4double energy = primary->GetParticleGun()->GetParticleEnergy();
  G4cout << "\n The run consists of " << NbOfEvents << " "<< Particle << " of "
         << G4BestUnit(energy,"Energy") << " through " 
	 << G4BestUnit(tickness,"Length") << " of "
	 << material->GetName() << " (density: " 
	 << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
  
  //frequency of processes
  G4int totalCount = 0;
  G4int survive = 0;  
  G4cout << "\n Process calls frequency --->";
  std::map<G4String,G4int>::iterator it;  
  for (it = procCounter.begin(); it != procCounter.end(); it++) {
     G4String procName = it->first;
     G4int    count    = it->second;
     totalCount += count; 
     G4cout << "\t" << procName << " = " << count;
     if (procName == "Transportation") survive = count;
  }  
  G4cout << G4endl;
  if (totalCount == 0) return;
  
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
  for (it = procCounter.begin(); it != procCounter.end(); it++) {
    G4String procName = it->first;  
    G4double massSigma = 
    emCalculator.GetCrossSectionPerVolume(energy,particle,
                                              procName,material)/density;
    if (particle == G4Gamma::Gamma())
       massSigma = 
       emCalculator.ComputeCrossSectionPerVolume(energy,particle,
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
	 	 
  //restore default format	 
  G4cout.precision(prec);         

  // remove all contents in procCounter 
  procCounter.clear();

  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
