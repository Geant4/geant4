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
// $Id: A01PhysicsList.cc,v 1.3 2002-12-13 11:34:34 gunter Exp $
// --------------------------------------------------------------
//


#include "A01PhysicsList.hh"

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include "g4std/iomanip"

#include "A01GeneralPhysics.hh"
#include "A01EMPhysics.hh"
#include "A01MuonPhysics.hh"
#include "A01HadronPhysics.hh"
#include "A01IonPhysics.hh"

A01PhysicsList::A01PhysicsList():  G4VModularPhysicsList()
{
  // default cut value  (1.0mm)
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);

  // General Physics
  RegisterPhysics( new A01GeneralPhysics("general") );

  // EM Physics
  RegisterPhysics( new A01EMPhysics("standard EM"));

  // Muon Physics
  RegisterPhysics(  new A01MuonPhysics("muon"));

   // Hadron Physics
  RegisterPhysics(  new A01HadronPhysics("hadron"));

  // Ion Physics
////////////////////////////////////////////  RegisterPhysics( new A01IonPhysics("ion"));


}

A01PhysicsList::~A01PhysicsList()
{
}

void A01PhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
  //   the default cut value for all particle types
  SetCutsWithDefault();
}



