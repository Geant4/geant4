//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4EmHighEnergyBuilder.cc,v 1.2 2004/11/24 13:18:02 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmHighEnergyBuilder
//
// Author:      V.Ivanchenko 03.05.2004
//
// Modified:
//
//----------------------------------------------------------------------------
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4EmHighEnergyBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4AnnihiToMuPair.hh"
#include "G4GammaConversionToMuons.hh"
#include "G4eeToHadrons.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmHighEnergyBuilder::G4EmHighEnergyBuilder(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmHighEnergyBuilder::~G4EmHighEnergyBuilder()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmHighEnergyBuilder::ConstructParticle()
{
  // Minimal set of particles
  G4Gamma::Gamma();
  G4Electron::Electron();
  G4Positron::Positron();
  G4MuonPlus::MuonPlus();
  G4MuonMinus::MuonMinus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmHighEnergyBuilder::ConstructProcess()
{
  // Add standard EM Processes for gamma
  G4ParticleDefinition* particle = G4Gamma::Gamma();
  G4ProcessManager* pmanager = particle->GetProcessManager();
  pmanager->AddDiscreteProcess( new G4GammaConversionToMuons() );

  // Add standard EM Processes for e+
  particle = G4Positron::Positron();
  pmanager = particle->GetProcessManager();
  G4AnnihiToMuPair* eplusmumu = new G4AnnihiToMuPair();
  pmanager->AddDiscreteProcess( eplusmumu );

  G4eeToHadrons* eplushadr = new G4eeToHadrons();
  pmanager->AddDiscreteProcess( eplushadr );
  
  // Built-in biasing factor for cross secsions
  eplusmumu->SetCrossSecFactor(1.0);
  eplushadr->SetCrossSecFactor(1.0);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

