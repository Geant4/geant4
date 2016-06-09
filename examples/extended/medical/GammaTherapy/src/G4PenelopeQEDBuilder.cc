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
//
// $Id: G4PenelopeQEDBuilder.cc,v 1.2 2006/06/29 17:27:46 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4PenelopeQEDBuilder
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

#include "G4PenelopeQEDBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4PenelopeCompton.hh"
#include "G4PenelopeGammaConversion.hh"
#include "G4PenelopePhotoElectric.hh"
#include "G4PenelopeRayleigh.hh"

#include "G4MultipleScattering.hh"

#include "G4PenelopeIonisation.hh"
#include "G4PenelopeBremsstrahlung.hh"
#include "G4PenelopeAnnihilation.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4PenelopeQEDBuilder::G4PenelopeQEDBuilder(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4PenelopeQEDBuilder::~G4PenelopeQEDBuilder()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4PenelopeQEDBuilder::ConstructParticle()
{
  G4Gamma::Gamma();
  G4Electron::Electron();
  G4Positron::Positron();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4PenelopeQEDBuilder::ConstructProcess()
{
  // Add standard EM Processes for gamma
  G4ParticleDefinition* particle = G4Gamma::Gamma();
  G4ProcessManager* pmanager = particle->GetProcessManager();

  pmanager->AddDiscreteProcess( new G4PenelopePhotoElectric() );
  pmanager->AddDiscreteProcess( new G4PenelopeCompton() );
  pmanager->AddDiscreteProcess( new G4PenelopeGammaConversion() );
  pmanager->AddDiscreteProcess( new G4PenelopeRayleigh() );

  // Add standard EM Processes for e-
  particle = G4Electron::Electron();
  pmanager = particle->GetProcessManager();

  pmanager->AddProcess(new G4MultipleScattering,     -1, 1,1);
  pmanager->AddProcess(new G4PenelopeIonisation,     -1, 2,2);
  pmanager->AddProcess(new G4PenelopeBremsstrahlung, -1,-1,3);

  // Add standard EM Processes for e+
  particle = G4Positron::Positron();
  pmanager = particle->GetProcessManager();

  pmanager->AddProcess(new G4MultipleScattering,    -1, 1,1);
  pmanager->AddProcess(new G4PenelopeIonisation,    -1, 2,2);
  pmanager->AddProcess(new G4PenelopeBremsstrahlung,-1,-1,3);
  pmanager->AddProcess(new G4PenelopeAnnihilation,   0,-1,4);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

