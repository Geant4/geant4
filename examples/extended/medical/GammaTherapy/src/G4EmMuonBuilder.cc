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
// $Id: G4EmMuonBuilder.cc,v 1.1 2004/12/02 10:34:08 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmMuonBuilder
//
// Author:      V.Ivanchenko 03.05.2004
//
// Modified:
// 24-11-2004 V.Ivanchenko Use the same radiation processes for mu+ and mu-
//
//----------------------------------------------------------------------------
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4EmMuonBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4MultipleScattering.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmMuonBuilder::G4EmMuonBuilder(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmMuonBuilder::~G4EmMuonBuilder()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmMuonBuilder::ConstructParticle()
{
  // Minimal set of particles
  G4Gamma::Gamma();
  G4Electron::Electron();
  G4Positron::Positron();
  G4MuonPlus::MuonPlus();
  G4MuonMinus::MuonMinus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmMuonBuilder::ConstructProcess()
{
  // Common processes for mu+ and mu-
  G4MultipleScattering* mumsc  = new G4MultipleScattering();
  G4MuIonisation*       muion  = new G4MuIonisation();
  G4MuBremsstrahlung*   mubrem = new G4MuBremsstrahlung();
  G4MuPairProduction*   mupair = new G4MuPairProduction();

  // Add standard EM Processes for mu+
  const G4ParticleDefinition* particle = G4MuonPlus::MuonPlus();
  G4ProcessManager* pmanager = particle->GetProcessManager();

  pmanager->AddProcess(mumsc,     -1, 1,1);
  pmanager->AddProcess(muion,     -1, 2,2);
  pmanager->AddProcess(mubrem,    -1,-1,3);
  pmanager->AddProcess(mupair,    -1,-1,4);


  // Add standard EM Processes for mu-
  particle = G4MuonMinus::MuonMinus();
  pmanager = particle->GetProcessManager();

  pmanager->AddProcess(mumsc,     -1, 1,1);
  pmanager->AddProcess(muion,     -1, 2,2);
  pmanager->AddProcess(mubrem,    -1,-1,3);
  pmanager->AddProcess(mupair,    -1,-1,4);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

