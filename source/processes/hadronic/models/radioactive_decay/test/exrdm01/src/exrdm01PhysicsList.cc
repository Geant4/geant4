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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		exrdm01PhysicsList.cc
//
// Version:		0.b.3
// Date:		29/02/00
// Author:		P R Truscott
// Organisation:	DERA UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		12115/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#include "exrdm01PhysicsList.hh"
#include "G4ParticleTypes.hh"
#include "G4IonConstructor.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4RadioactiveDecay.hh"
///////////////////////////////////////////////////////////////////////////////
//
exrdm01PhysicsList::exrdm01PhysicsList()
{;}
///////////////////////////////////////////////////////////////////////////////
//
exrdm01PhysicsList::~exrdm01PhysicsList()
{;}
///////////////////////////////////////////////////////////////////////////////
//
void exrdm01PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  G4IonTable *theIonTable =  (G4IonTable*) (G4ParticleTable::GetParticleTable()->GetIonTable());
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoE::NeutrinoEDefinition();
  G4Alpha::AlphaDefinition();

  G4Gamma::GammaDefinition();
  G4GenericIon::GenericIonDefinition();
}
///////////////////////////////////////////////////////////////////////////////
//
void exrdm01PhysicsList::ConstructProcess()
{
  //
  //
  // Only physical processes are transportation and radioactive decay.
  //
  AddTransportation();
  //
  //
  // Declare radioactive decay to the GenericIon in the IonTable.
  //
  G4IonTable *theIonTable =  (G4IonTable*) (G4ParticleTable::GetParticleTable()->
    GetIonTable());
  G4RadioactiveDecay *theRadioactiveDecay = new G4RadioactiveDecay();
  for (G4int i=0; i<theIonTable->Entries(); i++)
  {
    G4String particleName = theIonTable->GetParticle(i)->GetParticleName();
    if (particleName == "GenericIon")
    {
      G4ProcessManager* pmanager =
        theIonTable->GetParticle(i)->GetProcessManager();
      pmanager->SetVerboseLevel(0);
      pmanager ->AddProcess(theRadioactiveDecay);
      pmanager ->SetProcessOrdering(theRadioactiveDecay, idxPostStep);
      pmanager ->SetProcessOrdering(theRadioactiveDecay, idxAtRest);
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
//
void exrdm01PhysicsList::SetCuts()
{
  SetCutsWithDefault();
  SetCutValue (5.0*mm,"GenericIon");
}
///////////////////////////////////////////////////////////////////////////////
//
