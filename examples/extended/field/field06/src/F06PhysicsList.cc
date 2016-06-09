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
/// \file field/field06/src/F06PhysicsList.cc
/// \brief Implementation of the F06PhysicsList class
//
//
//
//

#include "F06PhysicsList.hh"

#include "F06ExtraPhysics.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4DecayPhysics.hh"
#include "G4ProcessTable.hh"

#include "G4PionDecayMakeSpin.hh"
#include "G4DecayWithSpin.hh"

#include "G4DecayTable.hh"
#include "G4MuonDecayChannelWithSpin.hh"
#include "G4MuonRadiativeDecayChannelWithSpin.hh"

F06PhysicsList::F06PhysicsList() : G4VModularPhysicsList() 
{
    RegisterPhysics(new G4DecayPhysics());
    RegisterPhysics(new F06ExtraPhysics());
}

F06PhysicsList::~F06PhysicsList()
{;}

void F06PhysicsList::ConstructParticle()
{
    G4Neutron::NeutronDefinition();
    G4Proton::ProtonDefinition();
    G4Electron::ElectronDefinition();
    G4AntiNeutrinoE::AntiNeutrinoEDefinition();
}

void F06PhysicsList::ConstructProcess()
{
    G4VModularPhysicsList::ConstructProcess();
}

void F06PhysicsList::SetCuts()
{
    SetCutsWithDefault();
}
