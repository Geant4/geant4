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
// $Id: G4EmHadronBuilder70.cc,v 1.1 2005-03-11 08:16:16 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmHadronBuilder70
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

#include "G4EmHadronBuilder70.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4MultipleScattering.hh"

#include "G4hIonisation70.hh"
#include "G4ionIonisation70.hh"

#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4GenericIon.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmHadronBuilder70::G4EmHadronBuilder70(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmHadronBuilder70::~G4EmHadronBuilder70()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmHadronBuilder70::ConstructParticle()
{
  G4Electron::Electron();
  G4Proton::Proton();
  G4PionPlus::PionPlus();
  G4PionMinus::PionMinus();
  G4KaonPlus::KaonPlus();
  G4KaonMinus::KaonMinus();
  G4GenericIon::GenericIon();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmHadronBuilder70::ConstructProcess()
{
  // Add standard EM Processes
  theParticleIterator->reset();
  G4hIonisation70* hion;
  const G4ParticleDefinition* pip = G4PionPlus::PionPlus();
  const G4ParticleDefinition* pin = G4PionMinus::PionMinus();

  while( (*theParticleIterator)() ){
    const G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    /*
    if(particle == pip) {
      pmanager->AddProcess(new G4MultipleScattering,-1,1,1); 
      hion = new G4hIonisation70();
      hion->SetBaseParticle(particle);
      pmanager->AddProcess(hion,-1,2,2);

    } else if(particle == pin) {
      pmanager->AddProcess(new G4MultipleScattering,-1,1,1); 
      hion = new G4hIonisation70();
      hion->SetBaseParticle(particle);
      pmanager->AddProcess(hion,-1,2,2);

    } if(particle == G4KaonPlus::KaonPlus()) {
      pmanager->AddProcess(new G4MultipleScattering,-1,1,1); 
      hion = new G4hIonisation70();
      hion->SetBaseParticle(pip);
      pmanager->AddProcess(hion,-1,2,2);

    } if(particle == G4KaonMinus::KaonMinus()) {
      pmanager->AddProcess(new G4MultipleScattering,-1,1,1); 
      hion = new G4hIonisation70();
      hion->SetBaseParticle(pin);
      pmanager->AddProcess(hion,-1,2,2);
  
    } else if(particle->GetPDGMass() > 110.*MeV) {
    */
    if(particle->GetPDGMass() > 110.*MeV) {
      G4String particleName = particle->GetParticleName();

      if (particleName == "GenericIon" || particleName == "alpha" || particleName == "He3") {

        pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);
        pmanager->AddProcess(new G4ionIonisation70,      -1, 2,2);

      } else if ((!particle->IsShortLived()) &&
	         (particle->GetPDGCharge() != 0.0) &&
	         (particle->GetParticleName() != "chargedgeantino")) {

        pmanager->AddProcess(new G4MultipleScattering,-1,1,1);
        pmanager->AddProcess(new G4hIonisation70,     -1,2,2);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

