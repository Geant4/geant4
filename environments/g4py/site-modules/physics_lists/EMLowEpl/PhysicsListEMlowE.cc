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
// $Id: PhysicsListEMlowE.cc,v 1.2 2010-06-21 12:29:42 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   PhysicsListEMlowE.cc
//
//                                         2006 Q
// ====================================================================
#include "PhysicsListEMlowE.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyRayleigh.hh"
#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4eMultipleScattering.hh"

// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////////////////////
PhysicsListEMlowE::PhysicsListEMlowE()
//////////////////////////////////////
{
  SetVerboseLevel(1);
  defaultCutValue = 1.*mm;   // default cut value  (1.0mm)
}


///////////////////////////////////////
PhysicsListEMlowE::~PhysicsListEMlowE()
///////////////////////////////////////
{
}


///////////////////////////////////////////
void PhysicsListEMlowE::ConstructParticle()
///////////////////////////////////////////
{
  G4Gamma::GammaDefinition();
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}


//////////////////////////////////////////
void PhysicsListEMlowE::ConstructProcess()
//////////////////////////////////////////
{
  AddTransportation();

  G4ProcessManager* pm;

  // ----------------------------------------------------------
  // gamma physics
  // ----------------------------------------------------------
  pm= G4Gamma::Gamma()-> GetProcessManager();
  pm-> AddDiscreteProcess(new G4LowEnergyPhotoElectric);
  pm-> AddDiscreteProcess(new G4LowEnergyCompton);
  pm-> AddDiscreteProcess(new G4LowEnergyGammaConversion);
  pm-> AddDiscreteProcess(new G4LowEnergyRayleigh);

  // ----------------------------------------------------------
  // electron physics
  // ----------------------------------------------------------
  G4LowEnergyIonisation* eion=       new G4LowEnergyIonisation;
  G4LowEnergyBremsstrahlung* ebrems= new G4LowEnergyBremsstrahlung;
  G4eMultipleScattering* msc=        new G4eMultipleScattering;

  pm= G4Electron::Electron()-> GetProcessManager();
  pm-> AddProcess(msc,    ordInActive,           1, 1);
  pm-> AddProcess(eion,   ordInActive,           2, 2);
  pm-> AddProcess(ebrems, ordInActive, ordInActive, 3);

  // ----------------------------------------------------------
  // positron physics
  // ----------------------------------------------------------
  eion=   new G4LowEnergyIonisation;
  ebrems= new G4LowEnergyBremsstrahlung;
  msc=    new G4eMultipleScattering;
  G4eplusAnnihilation* annihilation= new G4eplusAnnihilation;

  pm= G4Positron::Positron()-> GetProcessManager();
  pm-> AddProcess(msc,          ordInActive, 1,           1);
  pm-> AddProcess(eion,         ordInActive, 2,           2);
  pm-> AddProcess(ebrems,       ordInActive, ordInActive, 3);
  pm-> AddProcess(annihilation, 0,           ordInActive, 4);

}


/////////////////////////////////
void PhysicsListEMlowE::SetCuts()
/////////////////////////////////
{
  SetCutsWithDefault();
}

