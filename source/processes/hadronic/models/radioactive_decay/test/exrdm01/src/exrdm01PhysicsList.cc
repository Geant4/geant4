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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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
