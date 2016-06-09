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
// $Id: G4EmQEDBuilder52.cc,v 1.2 2005/05/08 16:12:53 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmQEDBuilder52
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

#include "G4EmQEDBuilder52.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4ComptonScattering52.hh"
#include "G4GammaConversion52.hh"
#include "G4PhotoElectricEffect52.hh"

#include "G4MultipleScattering52.hh"

#include "G4eIonisation52.hh"
#include "G4eBremsstrahlung52.hh"
#include "G4eplusAnnihilation52.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmQEDBuilder52::G4EmQEDBuilder52(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmQEDBuilder52::~G4EmQEDBuilder52()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmQEDBuilder52::ConstructParticle()
{
  G4Gamma::Gamma();
  G4Electron::Electron();
  G4Positron::Positron();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmQEDBuilder52::ConstructProcess()
{
  // Add standard EM Processes for gamma
  G4ParticleDefinition* particle = G4Gamma::Gamma();
  G4ProcessManager* pmanager = particle->GetProcessManager();

  pmanager->AddDiscreteProcess( new G4PhotoElectricEffect52() );
  pmanager->AddDiscreteProcess( new G4ComptonScattering52() );
  pmanager->AddDiscreteProcess( new G4GammaConversion52() );

  // Add standard EM Processes for e-
  particle = G4Electron::Electron();
  pmanager = particle->GetProcessManager();

  pmanager->AddProcess(new G4MultipleScattering52, -1, 1,1);
  pmanager->AddProcess(new G4eIonisation52,        -1, 2,2);
  pmanager->AddProcess(new G4eBremsstrahlung52,    -1,-1,3);

  // Add standard EM Processes for e+
  particle = G4Positron::Positron();
  pmanager = particle->GetProcessManager();

  pmanager->AddProcess(new G4MultipleScattering52, -1, 1,1);
  pmanager->AddProcess(new G4eIonisation52,        -1, 2,2);
  pmanager->AddProcess(new G4eBremsstrahlung52,    -1,-1,3);
  pmanager->AddProcess(new G4eplusAnnihilation52,   0,-1,4);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

