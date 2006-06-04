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
// $Id: PhysicsListEMstd.cc,v 1.3 2006-06-04 21:36:35 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   PhysicsListEMstd.cc
//
//                                         2006 Q
// ====================================================================
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
#include "G4MultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"


// ====================================================================
//
// class description
//
// ====================================================================

////////////////////////////////////
PhysicsListEMstd::PhysicsListEMstd()
////////////////////////////////////
{
  SetVerboseLevel(1);
  defaultCutValue = 1.*mm;   // default cut value  (1.0mm)
}


/////////////////////////////////////
PhysicsListEMstd::~PhysicsListEMstd()
/////////////////////////////////////
{
}


//////////////////////////////////////////
void PhysicsListEMstd::ConstructParticle()
//////////////////////////////////////////
{
  G4Gamma::GammaDefinition();
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}


/////////////////////////////////////////
void PhysicsListEMstd::ConstructProcess()
/////////////////////////////////////////
{
  AddTransportation();

  G4ProcessManager* pm;

  // ----------------------------------------------------------
  // gamma physics
  // ----------------------------------------------------------
  pm= G4Gamma::Gamma()-> GetProcessManager();
  pm-> AddDiscreteProcess(new G4PhotoElectricEffect);
  pm-> AddDiscreteProcess(new G4ComptonScattering);
  pm-> AddDiscreteProcess(new G4GammaConversion);

  // ----------------------------------------------------------
  // electron physics
  // ----------------------------------------------------------
  G4MultipleScattering* msc=    new G4MultipleScattering;
  G4eIonisation*        eion=   new G4eIonisation;
  G4eBremsstrahlung*    ebrems= new G4eBremsstrahlung;

  pm= G4Electron::Electron()->GetProcessManager();
  pm-> AddProcess(msc,    ordInActive,           1, 1);
  pm-> AddProcess(eion,   ordInActive,           2, 2);
  pm-> AddProcess(ebrems, ordInActive, ordInActive, 3);

  // ----------------------------------------------------------
  // positron physics
  // ----------------------------------------------------------
  msc=    new G4MultipleScattering;
  eion=   new G4eIonisation;
  ebrems= new G4eBremsstrahlung;
  G4eplusAnnihilation* annihilation= new G4eplusAnnihilation;

  pm= G4Positron::Positron()-> GetProcessManager();
  pm-> AddProcess(msc,          ordInActive, 1,           1);
  pm-> AddProcess(eion,         ordInActive, 2,           2);
  pm-> AddProcess(ebrems,       ordInActive, ordInActive, 3);
  pm-> AddProcess(annihilation, 0,           ordInActive, 4);

}


////////////////////////////////
void PhysicsListEMstd::SetCuts()
////////////////////////////////
{
  SetCutsWithDefault();
}

