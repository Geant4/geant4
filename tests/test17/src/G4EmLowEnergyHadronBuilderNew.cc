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
// $Id: G4EmLowEnergyHadronBuilderNew.cc,v 1.1 2004-05-26 11:39:09 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmLowEnergyHadronBuilderNew
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

#include "G4EmLowEnergyHadronBuilderNew.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4MultipleScattering.hh"

#include "G4hLowEnergyIonisationMA.hh"
#include "G4hLowEnergyIonisation.hh"
#include "G4EmProcessOptions.hh"

#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmLowEnergyHadronBuilderNew::G4EmLowEnergyHadronBuilderNew(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmLowEnergyHadronBuilderNew::~G4EmLowEnergyHadronBuilderNew()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmLowEnergyHadronBuilderNew::ConstructParticle()
{
  G4Electron::Electron();
  G4Proton::Proton();
  G4GenericIon::GenericIon();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmLowEnergyHadronBuilderNew::ConstructProcess()
{
  // Add standard EM Processes
  theParticleIterator->reset();

  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();


    if(particle->GetPDGMass() > 110.*MeV && !particle->IsShortLived() &&
       particle->GetPDGCharge() != 0.0) {

       pmanager->AddProcess(new G4MultipleScattering,-1,1,1);

       if(particle->GetPDGCharge() > 0.0) {
	 G4hLowEnergyIonisationMA* hion = new G4hLowEnergyIonisationMA();
         pmanager->AddProcess(hion,-1,2,2);
	 if(particle->GetPDGMass() > 900.*MeV) pmanager->AddProcess(&nstop,-1,3,-1);
       } else {
	 G4hLowEnergyIonisation* hion = new G4hLowEnergyIonisation();
         pmanager->AddProcess(hion,-1,2,2);
       }
    }
  }
  //G4EmProcessOptions emOptions;
  //emOptions.SetVerbose(2,"hLowEIoni");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

