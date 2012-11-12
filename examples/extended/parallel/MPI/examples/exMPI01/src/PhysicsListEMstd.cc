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
/// @file PhysicsListEMstd.cc
/// @brief Define standard EM physics list

#include "PhysicsListEMstd.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

// --------------------------------------------------------------------------
PhysicsListEMstd::PhysicsListEMstd()
  : G4VPhysicsConstructor("EM-std")
{
}

// --------------------------------------------------------------------------
PhysicsListEMstd::~PhysicsListEMstd()
{
}

// --------------------------------------------------------------------------
void PhysicsListEMstd::ConstructParticle()
{
}

// --------------------------------------------------------------------------
void PhysicsListEMstd::ConstructProcess()
{
  G4ProcessManager* pm;

  // ----------------------------------------------------------
  // gamma physics
  // ----------------------------------------------------------
  pm = G4Gamma::Gamma()-> GetProcessManager();
  pm-> AddDiscreteProcess(new G4PhotoElectricEffect);
  pm-> AddDiscreteProcess(new G4ComptonScattering);
  pm-> AddDiscreteProcess(new G4GammaConversion);

  // ----------------------------------------------------------
  // electron physics
  // ----------------------------------------------------------
  G4eMultipleScattering* msc =   new G4eMultipleScattering;
  G4eIonisation*        eion =   new G4eIonisation;
  G4eBremsstrahlung*    ebrems = new G4eBremsstrahlung;

  pm = G4Electron::Electron()->GetProcessManager();
  pm-> AddProcess(msc,    ordInActive,           1, 1);
  pm-> AddProcess(eion,   ordInActive,           2, 2);
  pm-> AddProcess(ebrems, ordInActive, ordInActive, 3);

  // ----------------------------------------------------------
  // positron physics
  // ----------------------------------------------------------
  msc =    new G4eMultipleScattering;
  eion =   new G4eIonisation;
  ebrems = new G4eBremsstrahlung;
  G4eplusAnnihilation* annihilation= new G4eplusAnnihilation;

  pm = G4Positron::Positron()-> GetProcessManager();
  pm-> AddProcess(msc,          ordInActive, 1,           1);
  pm-> AddProcess(eion,         ordInActive, 2,           2);
  pm-> AddProcess(ebrems,       ordInActive, ordInActive, 3);
  pm-> AddProcess(annihilation, 0,           ordInActive, 4);
}
