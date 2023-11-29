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
// S. Incerti (incerti@cenbg.in2p3.fr)
//

#include "G4EmDNAPhysics_option4.hh"
#include "G4EmDNABuilder.hh"
#include "G4SystemOfUnits.hh"

// ions
#include "G4Alpha.hh"
#include "G4DNAGenericIonsManager.hh"

// utilities
#include "G4EmParameters.hh"
#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"
#include "G4EmBuilder.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmDNAPhysics_option4);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDNAPhysics_option4::G4EmDNAPhysics_option4(G4int ver, const G4String& nam)
  : G4EmDNAPhysics(ver, nam)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysics_option4::ConstructProcess()
{
  // parameters
  G4EmParameters* param = G4EmParameters::Instance();
  const G4double emaxDNA = 1.*CLHEP::MeV;
  const G4double emaxIonDNA = 300.*CLHEP::MeV;
  const G4double eminBorn = 500.*CLHEP::keV;
  const G4bool fast = param->DNAFast();
  const G4bool st = param->DNAStationary();
  if(verboseLevel > 1) {
    G4cout << "### " << GetPhysicsName() 
	   << " Construct Processes EmaxDNA(MeV)= " 
           << emaxDNA/CLHEP::MeV << "; useMSC: " << fast 
           << "; stationary: " << st << G4endl;
  }
  G4DNAGenericIonsManager* genericIonsManager
    = G4DNAGenericIonsManager::Instance();

  // standard physics
  G4EmDNABuilder::ConstructStandardEmPhysics(emaxDNA, emaxIonDNA, 
                                             emaxIonDNA, emaxIonDNA,
                                             dnaGS, fast);

  // DNA physics
  G4EmDNABuilder::ConstructDNAElectronPhysics(emaxDNA, 4, fast, st);
  G4EmDNABuilder::ConstructDNAProtonPhysics(eminBorn, emaxIonDNA, 4, fast, st);
  G4EmDNABuilder::ConstructDNAIonPhysics(emaxIonDNA, st);

  G4ParticleDefinition* part = genericIonsManager->GetIon("hydrogen");
  G4EmDNABuilder::ConstructDNALightIonPhysics(part, 0, 4, emaxIonDNA, fast, st);

  part = G4Alpha::Alpha();
  G4EmDNABuilder::ConstructDNALightIonPhysics(part, 2, 4, emaxIonDNA, fast, st);

  part = genericIonsManager->GetIon("alpha+");
  G4EmDNABuilder::ConstructDNALightIonPhysics(part, 1, 4, emaxIonDNA, fast, st);

  part = genericIonsManager->GetIon("helium");
  G4EmDNABuilder::ConstructDNALightIonPhysics(part, 0, 4, emaxIonDNA, fast, st);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
