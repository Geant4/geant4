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
// $Id: G4EmLowEnergyHadronBuilder.cc,v 1.1 2004-05-26 11:39:09 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmLowEnergyHadronBuilder
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

#include "G4EmLowEnergyHadronBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4MultipleScattering.hh"

#include "G4hLowEnergyIonisation.hh"
#include "PhysicsListMessenger.hh"

#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmLowEnergyHadronBuilder::G4EmLowEnergyHadronBuilder(PhysicsListMessenger* messenger,
   const G4String& name)
   :  G4VPhysicsConstructor(name)
{
  if(messenger) messenger->SetLowEnergyHadronConstructor(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmLowEnergyHadronBuilder::~G4EmLowEnergyHadronBuilder()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmLowEnergyHadronBuilder::ConstructParticle()
{
  G4Electron::Electron();
  G4Proton::Proton();
  G4GenericIon::GenericIon();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmLowEnergyHadronBuilder::ConstructProcess()
{
  // Add standard EM Processes
  theParticleIterator->reset();

  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();


    if(particle->GetPDGMass() > 110.*MeV &&
       !particle->IsShortLived() &&
       particle->GetPDGCharge() != 0.0) {

	G4hLowEnergyIonisation* hion = new G4hLowEnergyIonisation();
	hio.push_back(hion);
        pmanager->AddProcess(new G4MultipleScattering,-1,1,1);
        pmanager->AddProcess(hion,-1,2,2);

    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmLowEnergyHadronBuilder::SetGammaCut(G4double cut)
{
  G4int n = hio.size();
  for(G4int i=0; i<n; i++) {
    (hio[i])->SetCutForSecondaryPhotons(cut);
    (hio[i])->SetFluorescence(true);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmLowEnergyHadronBuilder::SetAugerCut(G4double cut)
{
  G4int n = hio.size();
  for(G4int i=0; i<n; i++) {
    (hio[i])->SetCutForAugerElectrons(cut);
    (hio[i])->SetFluorescence(true);
    (hio[i])->ActivateAugerElectronProduction(true);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmLowEnergyHadronBuilder::ActivateFluorescence(G4bool val)
{
  G4int n = hio.size();
  for(G4int i=0; i<n; i++) {
    (hio[i])->SetFluorescence(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmLowEnergyHadronBuilder::ActivateAuger(G4bool val)
{
  G4int n = hio.size();
  for(G4int i=0; i<n; i++) {
    (hio[i])->ActivateAugerElectronProduction(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

